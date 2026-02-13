# Multi-Omics Integration Shiny App

A Shiny web application for integrating spatial proteomics data with single-cell RNA-seq data to infer disease-enriched or associated cell types within histopathological regions.

## Deployment

This app is deployed on AWS EC2 with Elastic IP: `98.90.240.106`

## Features

- **Step 1**: Upload spatial proteomics data and metadata
- **Step 2**: Identify region-specific proteins using limma differential expression
- **Step 3**: Upload scRNA-seq data (10X format)
- **Step 4**: Process and filter scRNA-seq data
- **Step 5**: Query method analysis for cell-type enrichment
- **Results**: Download comprehensive results and PDF reports

## Requirements

Install required R packages:

```r
install.packages(c(
  "shiny",
  "shinydashboard",
  "shinyjs",
  "DT",
  "dplyr",
  "readxl",
  "ggplot2",
  "reshape2",
  "RColorBrewer",
  "Matrix",
  "limma",
  "ggrepel",
  "VennDiagram",
  "grid"
))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Seurat")
```

## Running Locally

```r
shiny::runApp("app.R")
```

## File Upload Limits

The app supports large file uploads (up to 10GB) for scRNA-seq data matrices.

## References

- Wang et al. Proteomes, 2025, 13(1):3 — DOI: 10.3390/proteomes13010003
- Wang et al. Proteomes, 2025, 13(2):17 - DOI: 10.3390/proteomes13020017

## License

[Add your license here]
