library(shiny)
library(shinydashboard)
library(tidyverse)
library(DT)
library(ggplot2)
library(ggrepel)
library(tools)

# Load data function
load_deseq_data <- function() {
  DATA_DIR <- "data/"
  csv_files <- list.files(DATA_DIR, pattern = "*.csv", full.names = TRUE)
  result_list <- list()
  
  for (file_path in csv_files) {
    file_name <- file_path_sans_ext(basename(file_path))
    tryCatch({
      data <- read_csv(file_path, show_col_types = FALSE)
      required_cols <- c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
      if (all(required_cols %in% colnames(data))) {
        result_list[[file_name]] <- data
      }
    }, error = function(e) {
      warning(paste("Error reading", file_name, ":", e$message))
    })
  }
  
  return(result_list)
}

# Load metadata function
load_metadata_files <- function() {
  META_DIR <- "meta/"
  if (!dir.exists(META_DIR)) {
    return(list())
  }
  
  csv_files <- list.files(META_DIR, pattern = "*.csv", full.names = TRUE)
  result_list <- list()
  
  for (file_path in csv_files) {
    file_name <- file_path_sans_ext(basename(file_path))
    tryCatch({
      data <- read_csv(file_path, show_col_types = FALSE)
      result_list[[file_name]] <- list(
        data = data,
        path = file_path,
        size = file.size(file_path)
      )
    }, error = function(e) {
      warning(paste("Error reading metadata file", file_name, ":", e$message))
    })
  }
  
  return(result_list)
}

# Load raw data function
load_raw_data_files <- function() {
  RAW_DIR <- "raw/"
  if (!dir.exists(RAW_DIR)) {
    return(list())
  }
  
  csv_files <- list.files(RAW_DIR, pattern = "*.csv", full.names = TRUE)
  result_list <- list()
  
  for (file_path in csv_files) {
    file_name <- file_path_sans_ext(basename(file_path))
    tryCatch({
      data <- read_csv(file_path, show_col_types = FALSE)
      # Check if first column contains gene names/IDs
      if (ncol(data) > 1) {
        result_list[[file_name]] <- list(
          data = data,
          path = file_path,
          size = file.size(file_path),
          samples = colnames(data)[-1],  # Exclude first column (gene names)
          genes = data[[1]]  # First column should contain gene names
        )
      }
    }, error = function(e) {
      warning(paste("Error reading raw data file", file_name, ":", e$message))
    })
  }
  
  return(result_list)
}

# UI
ui <- dashboardPage(
  dashboardHeader(title = "ROS RNA-seq Dashboard"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Introduction", tabName = "intro", icon = icon("info-circle")),
      menuItem("Data Summary", tabName = "summary", icon = icon("table")),
      menuItem("Plots", tabName = "plots", icon = icon("chart-line")),
      menuItem("Gene Search", tabName = "gene_search", icon = icon("search")),
      menuItem("Boxplot Analysis", tabName = "boxplot", icon = icon("chart-bar")),
      menuItem("Metadata", tabName = "metadata", icon = icon("file-alt")),
      menuItem("Download", tabName = "download", icon = icon("download"))
    )
  ),
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "intro",
              fluidRow(
                box(
                  title = "Overview", status = "primary", solidHeader = TRUE, width = 12,  
                  
                  p("The project aims to unravel the mechanisms behind TKI resistance in ROS1+ lung cancer using patient-derived cell lines (PDCLs)."),
                  
                  strong("We'll focus on Transcriptomics:"),
                  tags$ul(
                    tags$li(
                      
                      " We study RNA expression to see which genes are turned on/off during TKI treatment.
                      Specifically, we'll look at changes in gene expression after treatment with crizotinib or entrectinib."
                    )
                  ),
                  
                  strong("The goal is to:"),
                  tags$ul(
                    tags$li("Understand which genes are up- or down-regulated upon TKI treatment."),
                    tags$li("Identify potential new resistance mechanisms or biomarkers for treatment response.")
                  ),
                  
                  br(),
                  p("The below columns must be included in the datasets:"),
                  tags$pre("gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj"),
                  br(),
                  textOutput("dataset_count")
                )
              )
      ),
      
      tabItem(tabName = "summary",
              fluidRow(
                column(
                  width = 12,
                  box(title = "Select Comparison", status = "primary", solidHeader = TRUE, width = 12,
                      uiOutput("comparison_selector")),
                  box(title = "Selected Table", status = "primary", solidHeader = TRUE, width = 12,
                      DT::dataTableOutput("selected_table"))
                )
              ),
              fluidRow(
                column(
                  width = 4,       # The width of your box (try 4, 6, or 8 as you like)
                  box(
                    title = "Summary Statistics", status = "primary", solidHeader = TRUE, width = 12,
                    tableOutput("summary_stats")
                  )
                )
              )
              
      ),
      
      tabItem(tabName = "plots",
              fluidRow(
                box(title = "Volcano Plot", status = "primary", solidHeader = TRUE, width = 12,
                    plotOutput("volcano_plot", height = "500px"))
              ),
              fluidRow(
                box(title = "MA Plot", status = "primary", solidHeader = TRUE, width = 12,
                    plotOutput("ma_plot", height = "500px"))
              )
      ),
      
      tabItem(tabName = "gene_search",
              fluidRow(
                box(
                  title = "Gene Search", status = "info", solidHeader = TRUE, width = 12,
                  fluidRow(
                    column(6, 
                           textInput("search_gene", "Enter Gene Name:", placeholder = "e.g., TP53, EGFR, ROS1")),
                    column(6,
                           br(),
                           actionButton("search_button", "Search", class = "btn-primary"))
                  ),
                  hr(),
                  conditionalPanel(
                    condition = "output.gene_found",
                    h4("Gene Expression Across Datasets"),
                    DT::dataTableOutput("gene_results_table")
                  ),
                  conditionalPanel(
                    condition = "!output.gene_found",
                    textOutput("gene_not_found_message")
                  )
                )
              ),
              fluidRow(
                conditionalPanel(
                  condition = "output.gene_found",
                  box(
                    title = "Gene Expression Comparison Plot", status = "info", solidHeader = TRUE, width = 12,
                    plotOutput("gene_comparison_plot", height = "400px")
                  )
                )
              )
      ),
      
      tabItem(tabName = "boxplot",
              fluidRow(
                box(
                  title = "Dynamic Boxplot Analysis", status = "success", solidHeader = TRUE, width = 12,
                  p("Create comparative boxplots from raw expression data by selecting control and treatment groups."),
                  
                  fluidRow(
                    column(4,
                           h4("1. Select Raw Data File"),
                           uiOutput("raw_file_selector"),
                           br(),
                           conditionalPanel(
                             condition = "output.raw_file_loaded",
                             h4("2. Select Control Samples"),
                             uiOutput("control_samples_selector")
                           )
                    ),
                    column(4,
                           conditionalPanel(
                             condition = "output.raw_file_loaded",
                             h4("3. Select Treatment Samples"),
                             uiOutput("treatment_samples_selector"),
                             br(),
                             h4("4. Select Gene"),
                             uiOutput("gene_selector_boxplot")
                           )
                    ),
                    column(4,
                           conditionalPanel(
                             condition = "output.boxplot_ready",
                             h4("5. Generate Plot"),
                             actionButton("generate_boxplot", "Create Boxplot", class = "btn-success"),
                             br(), br(),
                             div(id = "sample_summary",
                                 h5("Sample Summary:"),
                                 textOutput("sample_counts"))
                           )
                    )
                  )
                )
              ),
              fluidRow(
                conditionalPanel(
                  condition = "output.show_boxplot",
                  box(
                    title = "Expression Boxplot", status = "success", solidHeader = TRUE, width = 12,
                    plotOutput("dynamic_boxplot", height = "500px")
                  )
                )
              ),
              fluidRow(
                conditionalPanel(
                  condition = "output.show_boxplot",
                  box(
                    title = "Statistical Summary", status = "info", solidHeader = TRUE, width = 6,
                    tableOutput("boxplot_stats")
                  ),
                  box(
                    title = "Raw Data Table", status = "info", solidHeader = TRUE, width = 6,
                    DT::dataTableOutput("boxplot_data_table")
                  )
                )
              )
      ),
      
      tabItem(tabName = "metadata",
              fluidRow(
                box(
                  title = "Metadata Files", status = "warning", solidHeader = TRUE, width = 12,
                  p("Browse and download metadata CSV files from the meta/ directory."),
                  hr(),
                  conditionalPanel(
                    condition = "output.meta_files_available",
                    fluidRow(
                      column(6,
                             h4("Available Metadata Files:"),
                             uiOutput("meta_file_selector")),
                      column(6,
                             br(),
                             downloadButton("download_meta", "Download Selected File", class = "btn-warning"),
                             br(), br(),
                             textOutput("meta_file_info"))
                    ),
                    hr(),
                    h4("Preview of Selected File:"),
                    DT::dataTableOutput("meta_preview_table")
                  ),
                  conditionalPanel(
                    condition = "!output.meta_files_available",
                    div(class = "alert alert-info",
                        h4("No Metadata Files Found"),
                        p("No CSV files were found in the meta/ directory."))
                  )
                )
              )
      ),
      
      tabItem(tabName = "download",
              fluidRow(
                box(title = "Filter & Download", status = "success", solidHeader = TRUE, width = 12,
                    numericInput("padj_cutoff", "Adjusted p-value cutoff:", 0.05, min = 0, max = 1, step = 0.01),
                    numericInput("lfc_cutoff", "log2 Fold Change cutoff (absolute):", 1, min = 0, step = 0.1),
                    downloadButton("download_filtered", "Download Filtered Results (CSV)"),
                    br(),
                    textOutput("filter_summary"))
              )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  deseq_data <- reactive({
    load_deseq_data()
  })
  
  metadata_files <- reactive({
    load_metadata_files()
  })
  
  raw_data_files <- reactive({
    load_raw_data_files()
  })
  
  output$dataset_count <- renderText({
    paste("Loaded datasets:", length(deseq_data()), "files found")
  })
  
  output$comparison_selector <- renderUI({
    data_list <- deseq_data()
    if (length(data_list) > 0) {
      tagList(
        h4(paste("Found", length(data_list), "comparison files")),
        selectInput("selected_csv", "Choose comparison:", 
                    choices = names(data_list), 
                    selected = names(data_list)[1])
      )
    } else {
      tagList(
        h4("No data found"),
        p("No CSV files found in data/ directory.")
      )
    }
  })
  
  output$selected_table <- renderDataTable({
    req(input$selected_csv)
    data_list <- deseq_data()
    df <- data_list[[input$selected_csv]]
    datatable(df, options = list(pageLength = 10, scrollX = TRUE,scrollY = "350px"))
  })
  
  output$summary_stats <- renderTable({
    req(input$selected_csv)
    df <- deseq_data()[[input$selected_csv]]
    tibble(
      Metric = c("Up-regulated genes (padj<0.05, log2FC>1)", 
                 "Down-regulated genes (padj<0.05, log2FC<-1)",
                 "Not significant",
                 "NA values (filtered out)",
                 "Total genes analyzed"),
      Count = c(
        sum(df$padj < 0.05 & df$log2FoldChange > 1, na.rm = TRUE),
        sum(df$padj < 0.05 & df$log2FoldChange < -1, na.rm = TRUE),
        sum(df$padj >= 0.05 | abs(df$log2FoldChange) <= 1, na.rm = TRUE),
        sum(is.na(df$padj)),
        nrow(df)
      )
    )
  })
  
  output$volcano_plot <- renderPlot({
    req(input$selected_csv)
    df <- deseq_data()[[input$selected_csv]] %>%
      drop_na(padj) %>%
      mutate(sig = padj < 0.05 & abs(log2FoldChange) > 1,
             regulation = case_when(
               padj < 0.05 & log2FoldChange > 1 ~ "Up-regulated",
               padj < 0.05 & log2FoldChange < -1 ~ "Down-regulated",
               TRUE ~ "Not significant"
             ))
    
    # Select genes with highest absolute log fold change values (that are also significant)
    label_genes <- df %>% 
      filter(sig) %>% 
      arrange(desc(abs(log2FoldChange))) %>% 
      head(10)
    
    ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
      geom_point(alpha = 0.6) +
      geom_text_repel(data = label_genes, aes(label = gene), box.padding = 0.5, max.overlaps = 15) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgray") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgray") +
      scale_color_manual(values = c("Not significant" = "grey", "Up-regulated" = "red", "Down-regulated" = "blue")) +
      theme_minimal() +
      labs(title = paste("Volcano plot:", input$selected_csv), 
           x = "log2 Fold Change", y = "-log10(padj)", color = "Regulation Status")
  })
  
  output$ma_plot <- renderPlot({
    req(input$selected_csv)
    df <- deseq_data()[[input$selected_csv]] %>%
      drop_na(padj) %>%
      mutate(sig = padj < 0.05 & abs(log2FoldChange) > 1)
    
    # Select genes with highest absolute log fold change values (that are also significant)
    label_genes <- df %>% 
      filter(sig) %>% 
      arrange(desc(abs(log2FoldChange))) %>% 
      head(10)
    
    # baseMean is the average of the
    # normalized counts for each gene, calculated across all samples in the dataset, 
    # regardless of group assignment
    # a general indicator of how highly a gene is expressed overall
    ggplot(df, aes(x = log10(baseMean), y = log2FoldChange, color = sig)) +
      geom_point(alpha = 0.6) +
      geom_text_repel(data = label_genes, aes(label = gene), box.padding = 0.5, max.overlaps = 15) +
      geom_hline(yintercept = 0, color = "darkgray") +
      geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "darkgray") +
      scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "purple")) +
      theme_minimal() +
      labs(title = paste("MA plot:", input$selected_csv), 
           x = "log10(Mean Expression)", y = "log2 Fold Change")
  })
  
  # Gene search functionality
  gene_search_results <- eventReactive(input$search_button, {
    req(input$search_gene)
    
    gene_name <- toupper(trimws(input$search_gene))  # Convert to uppercase and trim whitespace
    data_list <- deseq_data()
    
    if (length(data_list) == 0) {
      return(NULL)
    }
    
    results <- list()
    
    for (dataset_name in names(data_list)) {
      df <- data_list[[dataset_name]]
      # Search for gene (case-insensitive)
      gene_row <- df[toupper(df$gene) == gene_name, ]
      
      if (nrow(gene_row) > 0) {
        gene_row$dataset <- dataset_name
        results[[dataset_name]] <- gene_row
      }
    }
    
    if (length(results) > 0) {
      # Combine all results
      combined_results <- bind_rows(results) %>%
        select(dataset, gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
        arrange(padj)
      
      return(combined_results)
    } else {
      return(data.frame())  # Return empty data frame if gene not found
    }
  })
  
  output$gene_found <- reactive({
    results <- gene_search_results()
    return(!is.null(results) && nrow(results) > 0)
  })
  outputOptions(output, "gene_found", suspendWhenHidden = FALSE)
  
  output$gene_results_table <- renderDataTable({
    results <- gene_search_results()
    req(nrow(results) > 0)
    
    datatable(results, 
              options = list(pageLength = 10, scrollX = TRUE),
              rownames = FALSE) %>%
      formatRound(columns = c("baseMean", "log2FoldChange", "lfcSE", "stat"), digits = 3) %>%
      formatSignif(columns = c("pvalue", "padj"), digits = 3)
  })
  
  output$gene_not_found_message <- renderText({
    results <- gene_search_results()
    if (!is.null(results) && nrow(results) == 0) {
      paste("Gene '", input$search_gene, "' not found in any dataset.")
    } else if (is.null(results)) {
      "No datasets loaded or search not performed."
    } else {
      ""
    }
  })
  
  output$gene_comparison_plot <- renderPlot({
    results <- gene_search_results()
    req(nrow(results) > 0)
    
    # Create significance categories
    plot_data <- results %>%
      mutate(
        significance = case_when(
          is.na(padj) ~ "Not available",
          padj < 0.001 ~ "p < 0.001",
          padj < 0.01 ~ "p < 0.01",
          padj < 0.05 ~ "p < 0.05",
          TRUE ~ "Not significant"
        ),
        regulation = case_when(
          is.na(padj) ~ "Not available",
          padj < 0.05 & log2FoldChange > 0 ~ "Up-regulated",
          padj < 0.05 & log2FoldChange < 0 ~ "Down-regulated",
          TRUE ~ "Not significant"
        )
      )
    
    # Create the plot
    p1 <- ggplot(plot_data, aes(x = reorder(dataset, log2FoldChange), y = log2FoldChange, fill = regulation)) +
      geom_col() +
      geom_text(aes(label = paste0("padj: ", signif(padj, 3))), 
                hjust = ifelse(plot_data$log2FoldChange >= 0, -0.1, 1.1),
                size = 3) +
      scale_fill_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue", 
                                   "Not significant" = "grey", "Not available" = "lightgrey")) +
      coord_flip() +
      theme_minimal() +
      labs(
        title = paste("Expression of", toupper(input$search_gene), "across datasets"),
        x = "Dataset",
        y = "log2 Fold Change",
        fill = "Regulation Status"
      ) +
      theme(axis.text.y = element_text(size = 10))
    
    print(p1)
  })
  
  # Metadata functionality
  output$meta_files_available <- reactive({
    meta_data <- metadata_files()
    return(length(meta_data) > 0)
  })
  outputOptions(output, "meta_files_available", suspendWhenHidden = FALSE)
  
  output$meta_file_selector <- renderUI({
    meta_data <- metadata_files()
    if (length(meta_data) > 0) {
      selectInput("selected_meta", "Choose metadata file:", 
                  choices = names(meta_data), 
                  selected = names(meta_data)[1])
    }
  })
  
  output$meta_file_info <- renderText({
    req(input$selected_meta)
    meta_data <- metadata_files()
    selected_meta <- meta_data[[input$selected_meta]]
    
    file_size_kb <- round(selected_meta$size / 1024, 2)
    n_rows <- nrow(selected_meta$data)
    n_cols <- ncol(selected_meta$data)
    
    paste0("File: ", input$selected_meta, ".csv\n",
           "Size: ", file_size_kb, " KB\n",
           "Dimensions: ", n_rows, " rows × ", n_cols, " columns")
  })
  
  output$meta_preview_table <- renderDataTable({
    req(input$selected_meta)
    meta_data <- metadata_files()
    selected_data <- meta_data[[input$selected_meta]]$data
    
    datatable(selected_data, 
              options = list(pageLength = 10, scrollX = TRUE, scrollY = "300px"),
              rownames = FALSE)
  })
  
  output$download_meta <- downloadHandler(
    filename = function() {
      paste0(input$selected_meta, ".csv")
    },
    content = function(file) {
      meta_data <- metadata_files()
      selected_data <- meta_data[[input$selected_meta]]$data
      write.csv(selected_data, file, row.names = FALSE)
    }
  )
  
  # Boxplot functionality
  output$raw_file_selector <- renderUI({
    raw_data <- raw_data_files()
    if (length(raw_data) > 0) {
      selectInput("selected_raw_file", "Choose raw data file:", 
                  choices = names(raw_data), 
                  selected = names(raw_data)[1])
    } else {
      div(class = "alert alert-warning",
          p("No CSV files found in raw/ directory."))
    }
  })
  
  output$raw_file_loaded <- reactive({
    raw_data <- raw_data_files()
    return(length(raw_data) > 0 && !is.null(input$selected_raw_file))
  })
  outputOptions(output, "raw_file_loaded", suspendWhenHidden = FALSE)
  
  output$control_samples_selector <- renderUI({
    req(input$selected_raw_file)
    raw_data <- raw_data_files()
    samples <- raw_data[[input$selected_raw_file]]$samples
    
    checkboxGroupInput("control_samples", 
                       "Select control samples:",
                       choices = samples,
                       selected = NULL)
  })
  
  output$treatment_samples_selector <- renderUI({
    req(input$selected_raw_file)
    raw_data <- raw_data_files()
    samples <- raw_data[[input$selected_raw_file]]$samples
    
    checkboxGroupInput("treatment_samples", 
                       "Select treatment samples:",
                       choices = samples,
                       selected = NULL)
  })
  
  
  output$gene_selector_boxplot <- renderUI({
    req(input$selected_raw_file)
    raw_data <- raw_data_files()
    genes <- raw_data[[input$selected_raw_file]]$genes
    
    selectizeInput("selected_gene_boxplot", 
                   "Search and select gene:",
                   choices = head(genes, 100),  # Only first 100 genes
                   selected = NULL,
                   options = list(
                     placeholder = "Type gene name to search...",
                     maxOptions = 100,
                     create = TRUE,  # Allow user to type gene names not in the list
                     selectOnTab = TRUE
                   ))
  })
  

  
  output$boxplot_ready <- reactive({
    return(!is.null(input$control_samples) && 
             !is.null(input$treatment_samples) && 
             !is.null(input$selected_gene_boxplot) &&
             length(input$control_samples) > 0 &&
             length(input$treatment_samples) > 0 &&
             input$selected_gene_boxplot != "")
  })
  outputOptions(output, "boxplot_ready", suspendWhenHidden = FALSE)
  
  output$sample_counts <- renderText({
    req(input$control_samples, input$treatment_samples)
    paste0("Control samples: ", length(input$control_samples), "\n",
           "Treatment samples: ", length(input$treatment_samples))
  })
  
  # Reactive values for boxplot
  boxplot_data <- eventReactive(input$generate_boxplot, {
    req(input$selected_raw_file, input$control_samples, input$treatment_samples, input$selected_gene_boxplot)
    
    raw_data <- raw_data_files()
    selected_data <- raw_data[[input$selected_raw_file]]$data
    
    # Find the gene row
    gene_row_idx <- which(selected_data[[1]] == input$selected_gene_boxplot)
    
    if (length(gene_row_idx) == 0) {
      return(NULL)
    }
    
    gene_data <- selected_data[gene_row_idx, ]
    
    # Prepare data for plotting
    all_samples <- c(input$control_samples, input$treatment_samples)
    
    plot_data <- data.frame(
      sample = all_samples,
      expression = as.numeric(gene_data[1, all_samples]),
      group = c(rep("Control", length(input$control_samples)),
                rep("Treatment", length(input$treatment_samples))),
      stringsAsFactors = FALSE
    )
    
    # Remove any samples with NA values
    plot_data <- plot_data[!is.na(plot_data$expression), ]
    
    return(list(
      plot_data = plot_data,
      gene_name = input$selected_gene_boxplot
    ))
  })
  
  output$show_boxplot <- reactive({
    data <- boxplot_data()
    return(!is.null(data) && nrow(data$plot_data) > 0)
  })
  outputOptions(output, "show_boxplot", suspendWhenHidden = FALSE)
  
  output$dynamic_boxplot <- renderPlot({
    data <- boxplot_data()
    req(data)
    
    plot_data <- data$plot_data
    gene_name <- data$gene_name
    
    # Calculate statistics for annotation
    control_data <- plot_data[plot_data$group == "Control", "expression"]
    treatment_data <- plot_data[plot_data$group == "Treatment", "expression"]
    
    # Perform t-test if both groups have enough samples
    # The p-value tells you if the difference is statistically significant.
    #Cohen's d tells you how meaningful or large that difference is in practical terms
    
    p_value <- NULL
    cohens_d <- NULL
    if (length(control_data) >= 2 && length(treatment_data) >= 2) {
      tryCatch({
        t_test_result <- t.test(treatment_data, control_data)
        p_value <- t_test_result$p.value
        
        # Calculate Cohen's d
        mean_diff <- mean(treatment_data) - mean(control_data)
        pooled_sd <- sqrt(((sd(treatment_data)^2 + sd(control_data)^2) / 2))
        cohens_d <- mean_diff / pooled_sd
      }, error = function(e) {
        p_value <- NULL
        cohens_d <- NULL
      })
    }
    
    # Create the plot
    p <- ggplot(plot_data, aes(x = group, y = expression, fill = group)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
      scale_fill_manual(values = c("Control" = "lightblue", "Treatment" = "lightcoral")) +
      theme_minimal() +
      labs(
        title = paste("Expression of", gene_name),
        subtitle = if (!is.null(p_value) && !is.null(cohens_d)) {
          paste0("p-value: ", signif(p_value, 3), 
                 ", Cohen's d: ", signif(cohens_d, 3))
        } else if (!is.null(p_value)) {
          paste("p-value:", signif(p_value, 3))
        } else "",
        x = "Group",
        y = "Expression Level",
        fill = "Group"
      ) +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none"
      )
    
    # Add significance annotation if p-value is available
    if (!is.null(p_value)) {
      max_y <- max(plot_data$expression, na.rm = TRUE)
      y_position <- max_y * 1.1
      
      significance <- ifelse(p_value < 0.001, "***",
                             ifelse(p_value < 0.01, "**",
                                    ifelse(p_value < 0.05, "*", "not significant")))
      
      p <- p + 
        annotate("text", x = 1.5, y = y_position, 
                 label = significance, size = 6) +
        annotate("segment", x = 1, xend = 2, y = y_position * 0.95, yend = y_position * 0.95)
    }
    
    print(p)
    
  })
  
  output$boxplot_stats <- renderTable({
    data <- boxplot_data()
    req(data)
    
    plot_data <- data$plot_data
    
    stats <- plot_data %>%
      group_by(group) %>%
      summarise(
        n = n(),
        mean = round(mean(expression, na.rm = TRUE), 3),
        median = round(median(expression, na.rm = TRUE), 3),
        sd = round(sd(expression, na.rm = TRUE), 3),
        min = round(min(expression, na.rm = TRUE), 3),
        max = round(max(expression, na.rm = TRUE), 3),
        .groups = 'drop'
      )
    
    return(stats)
  })
  
  output$boxplot_data_table <- renderDataTable({
    data <- boxplot_data()
    req(data)
    
    plot_data <- data$plot_data %>%
      arrange(group, sample)
    
    datatable(plot_data, 
              options = list(pageLength = 10, scrollX = TRUE),
              rownames = FALSE) %>%
      formatRound(columns = "expression", digits = 3)
  })
  
  output$download_filtered <- downloadHandler(
    filename = function() {
      paste0("filtered_", input$selected_csv, "_padj", input$padj_cutoff, "_FC", input$lfc_cutoff, ".csv")
    },
    content = function(file) {
      df <- deseq_data()[[input$selected_csv]] %>%
        filter(!is.na(padj), padj < input$padj_cutoff, abs(log2FoldChange) > input$lfc_cutoff) %>%
        arrange(padj)
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  output$filter_summary <- renderText({
    req(input$selected_csv)
    df <- deseq_data()[[input$selected_csv]]
    filtered_count <- sum(!is.na(df$padj) & df$padj < input$padj_cutoff & abs(df$log2FoldChange) > input$lfc_cutoff)
    paste0("Current filter: padj < ", input$padj_cutoff,
           " and |log2FC| > ", input$lfc_cutoff,
           "\nNumber of genes passing filter: ", filtered_count,
           " out of ", nrow(df), " total genes.")
  })
}

# Run the app
shinyApp(ui = ui, server = server)