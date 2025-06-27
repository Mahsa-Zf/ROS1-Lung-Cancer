# ROS1+ Lung Cancer

ROS1 is a gene that encodes a receptor tyrosine kinase. Under normal circumstances, this protein helps regulate cell growth and division. However, in ROS1+ lung cancer, the ROS1 gene becomes abnormally fused with another gene (a fusion partner), leading to uncontrolled cell growth.

This repository contains a jupyter notebook with an R kernel for analysis and an interactive Shiny dashboard for the visualization of the results of the RNA-seq data analysis from patient-derived cell lines (PDCLs). The goal is to investigate mechanisms of tyrosine kinase inhibitor (TKI) resistance in ROS1-positive lung cancer by assessing the up/down regulated genes.

## Project Overview

The ros Jupyter Notebook using an **R kernel** explores, pre-processes and prepares the data related to ROS1+ lung cancer from 2 studies. 

The studies used for the analysis are:
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239844
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214715

**The raw count data were directly downloaded from these pages.**

The dashboard is designed to use the prepared data to interactively:

- Explore transcriptomic changes in ROS1+ lung cancer cell lines after treatment with TKIs (crizotinib or entrectinib).
- Identify genes that are up- or down-regulated in response to TKI treatment.
- Facilitate discovery of new resistance mechanisms or potential biomarkers for treatment response.


## Features

- **Data Summary:** View and interact with summary tables of differential gene expression results.
- **Plots:** Generate volcano and MA plots for visualizing significant gene expression changes.
- **Gene Search:** Search for specific genes across all loaded datasets and compare their expression patterns.
- **Boxplot Analysis:** Create comparative boxplots for selected genes between control and treatment groups using raw expression data.
- **Metadata Browser:** Browse and download associated metadata files.
- **Download Module:** Filter and download results based on significance thresholds.

**Data Summary Columns Meaning:**
- BaseMean: average normalized count value for this gene across all samples. It indicates the gene's overall expression level
- log2FoldChange: This represents the effect size or how much expression changed with treatment. A value of 1 means the gene is 2 times more highly expressed in your treatment compared to the reference i.e. control.
- IfcSE: standard error of the log2 fold change estimate.
- stat: Wald statistic for the hypothesis test. Higher values indicate stronger evidence against the null hypothesis.
- padj: adjusted p-value for multiple testing. It indicates the significance of the result after correcting for multiple comparisons.

## Run the app:

https://mzamani.shinyapps.io/ROS-Seq/

## Contact

For questions or contributions, please open an issue or contact the repository maintainer.

**Note:**
The file `Gene expression in ROS1.docx` provides additional background information and is not directly used by the dashboard application.

