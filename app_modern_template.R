options(shiny.maxRequestSize = 10000 * 1024^2)


# Multi-Omics Integration Shiny App
# Integrating Spatial Proteomics and scRNA-seq
# Based on: https://github.com/outobi/query_and_overlap_method_pipeline

library(shiny)
library(bslib)          
library(bsicons)        
library(DT)
library(dplyr)
library(readxl)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(Matrix)
library(Seurat)
library(limma)
library(ggrepel)
library(VennDiagram)
library(grid)

# 1. 定义低明度、高质感的深蓝主题
my_theme <- bs_theme(
  version = 5,
  # 官方 AbbVie 深蓝色调
  primary = "#002F6C",
  # 辅助色改为中性灰，防止干扰
  secondary = "#6C757D",
  "navbar-bg" = "#ffffff",
  "navbar-color" = "#000000"
)


# UI Definition
ui <- page_navbar(
  theme = my_theme,
  title = "Multi-Omics Integration Pipeline",
  id = "navbar",
  fillable = FALSE,
  header = tags$head(
    tags$style(HTML("
      .bg-primary,
      .text-bg-primary,
      .btn-primary {
        background-color: #002F6C !important;
        border-color: #002F6C !important;
      }
      .btn-primary:hover,
      .btn-primary:focus {
        background-color: #001f4a !important;
        border-color: #001f4a !important;
      }
    "))
  ),

  # --- Home Tab ---
  nav_panel(
    title = "Home", 
    value = "home",
    icon = bs_icon("house-fill"),
    layout_column_wrap(
      width = 1,
      card(
        card_header(class = "bg-primary text-white", "🧬 Spatial Proteomics Integration with scRNA-seq transcriptomics"),
        div(
          h3("Welcome to Multi-Omics Integration Pipeline"),
          p("This pipeline integrates spatial proteomics data (laser microdissection LC/MS proteomics, Nanostring GeoMx DSP) with single-cell RNA-seq data 
          to infer disease-enriched or associated cell types within histopathological regions."),
          hr(),
          h4("📊 Analysis Workflow:"),
          tags$ol(
            tags$li(tags$b("Proteomics Upload:"), " Upload your spatial proteomics expression matrix and sample metadata"),
            tags$li(tags$b("Region-Specific Proteins:"), " Identify upregulated region-specific proteins using statistical tests"),
            tags$li(tags$b("scRNA-seq Upload:"), " Upload your single-cell RNA-seq data"),
            tags$li(tags$b("scRNA-seq Processing:"), " Process and validate scRNA-seq data with cell type annotations"),
            tags$li(tags$b("Query Method:"), " Compute cell-type specific enrichment z-scores of region-specific proteins"),
            tags$li(tags$b("Results:"), " Download results, plots, and summary tables")
          ),
          hr(),
          div(
            style = "background-color: #ffffff; padding: 15px; border: 1px solid #dee2e6; border-radius: 5px;",
            h4("📚 References:"),
            p("• Wang et al. Proteomes, 2025, 13(1):3 — DOI: 10.3390/proteomes13010003"),
            p("• Wang et al. Proteomes, 2025, 13(2):17 - DOI: 10.3390/proteomes13020017")
          ),
          hr(),
          actionButton("start_btn", "Get Started →", class = "btn btn-primary btn-lg")
        )
      )
    )
  ),

  # --- Step 1: Proteomics Upload ---
  nav_panel(
    title = "Step 1: Proteomics", 
    value = "proteomics_upload",
    icon = bs_icon("upload"),
    layout_column_wrap(
      width = 1/2,
      card(
        card_header(class = "bg-primary text-white", "1. Spatial Proteomics Expression Matrix"),
        p("Upload your spatial proteomics expression matrix (proteins × samples)."),
        div(
          class = "alert alert-light border",
          p(tags$b("Expected Format:"), style = "margin-bottom: 5px;"),
          tags$table(
            class = "table table-sm table-bordered",
            style = "font-family: monospace; font-size: 11px;",
            tags$tr(
              tags$th("Protein"),
              tags$th("S1"),
              tags$th("S2"),
              tags$th("S3"),
              tags$th("...")
            ),
            tags$tr(
              tags$td("CD3E"),
              tags$td("12.206"),
              tags$td("13.397"),
              tags$td("12.384"),
              tags$td("...")
            ),
            tags$tr(
              tags$td("VIM"),
              tags$td("13.963"),
              tags$td("12.083"),
              tags$td("13.189"),
              tags$td("...")
            ),
            tags$tr(
              tags$td("..."),
              tags$td("..."),
              tags$td("..."),
              tags$td("..."),
              tags$td("...")
            )
          ),
          p(class = "text-muted small", style = "margin-top: 10px;",
            "• First column: Protein names/IDs (gene symbols)", tags$br(),
            "• Remaining columns: Sample IDs (matching metadata)", tags$br(),
            "• Values: Protein abundance (log2 transformed, no missing values)")
        ),
        div(
          class = "alert alert-warning",
          p(bs_icon("lightbulb"), tags$b(" Test with example data:")),
          p("Laser microdissection coupled LC/MS spatial proteomics expression matrix"),
          p("Stroma region (n = 11) vs Tumor region (n = 12) in human TNBC FFPE tissue sections"),
          p("6842 protein IDs"),
          downloadButton("download_example_expr", "Download Example Data", class = "btn-warning btn-sm")
        ),
        fileInput("spatial_expr_file", "Upload Expression Matrix:", accept = c(".csv", ".xlsx", ".tsv")),
        uiOutput("spatial_expr_preview")
      ),
      card(
        card_header(class = "bg-primary text-white", "2. Sample Metadata"),
        p("Upload sample metadata with region/group annotations."),
        div(
          class = "alert alert-light border",
          p(tags$b("Expected Format:"), style = "margin-bottom: 5px;"),
          tags$table(
            class = "table table-sm table-bordered",
            style = "font-family: monospace; font-size: 11px;",
            tags$tr(
              tags$th("Sample_ID"),
              tags$th("Sample_group"),
              tags$th("Other_Info")
            ),
            tags$tr(
              tags$td("S1"),
              tags$td("Tumor"),
              tags$td("...")
            ),
            tags$tr(
              tags$td("S2"),
              tags$td("Stroma"),
              tags$td("...")
            ),
            tags$tr(
              tags$td("..."),
              tags$td("..."),
              tags$td("...")
            )
          ),
          p(class = "text-muted small", style = "margin-top: 10px;",
            "• First column: Sample ID (matching expression matrix)", tags$br(),
            "• Second column: Sample group for region-specific differential analysis", tags$br(),
            "• Additional columns: Optional (e.g., Patient_ID, Tissue_Type, Clinical_Info)"
          )
        ),
        div(
          class = "alert alert-warning",
          downloadButton("download_example_metadata", "Download Example Metadata", class = "btn-warning btn-sm")
        ),
        fileInput("metadata_file", "Upload Metadata:", accept = c(".csv", ".xlsx", ".tsv")),
        uiOutput("metadata_preview")
      )
    ),
    card(
      actionButton("validate_proteomics_btn", "Validate Proteomics Data", class = "btn btn-primary btn-lg"),
      hr(),
      verbatimTextOutput("proteomics_validation_summary")
    )
  ),

  # --- Step 2: Region Analysis ---
  nav_panel(
    title = "Step 2: Regions", 
    value = "region_genes",
    icon = bs_icon("search"),
    card(
      card_header(class = "bg-primary text-white", "Region-Specific Protein Identification"),
      div(
        class = "alert alert-info",
        p(bs_icon("info-circle"), " This step identifies proteins significantly upregulated in each region using limma differential expression analysis.")
      )
    ),
    layout_sidebar(
      sidebar = sidebar(
        title = "Analysis Parameters",
        selectInput("group_column", "Region Column:", choices = NULL),
        numericInput("pval_threshold", "Adj. P-Value Threshold:", value = 0.05, min = 0, max = 1, step = 0.01),
        numericInput("fc_threshold", "Fold Change Threshold:", value = 2, min = 1, max = 10, step = 0.5),
        actionButton("run_region_analysis_btn", "Run Analysis", class = "btn btn-primary w-100")
      ),
      navset_card_pill(
        id = "region_tabs",
        nav_panel(
          "PCA Plot",
          plotOutput("pca_plot", height = "600px"),
          downloadButton("download_pca_plot", "Download Plot", class = "btn-sm")
        ),
        nav_panel(
          "Venn Diagram",
          plotOutput("venn_plot", height = "600px"),
          verbatimTextOutput("overlap_stats"),
          downloadButton("download_venn_plot", "Download Plot", class = "btn-sm")
        ),
        nav_panel(
          "DE Results Table",
          DTOutput("region_genes_table"),
          downloadButton("download_region_genes", "Download Table", class = "btn-sm")
        )
      )
    ),
    uiOutput("region_analysis_status")
  ),

  # --- Step 3: scRNA-seq Upload ---
  nav_panel(
    title = "Step 3: scRNA-seq", 
    value = "scrna_upload",
    icon = bs_icon("upload"),
    layout_column_wrap(
      width = 1/2,
      # Column 1: Matrix & Genes
      card(
        card_header(class = "bg-primary text-white", "1. Sparse Matrix Expression File"),
        p("Upload the count matrix in sparse Matrix Market format (.mtx)."),
        
        div(style = "padding: 15px; background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 5px; margin-bottom: 15px;",
          tags$strong("Expected Format:"),
         
          tags$ul(style = "margin-top: 10px; font-size: 0.9em; color: #6c757d;",
            
            tags$li("Rows = genes, Columns = cells"),
            tags$li("Can be gzipped (.mtx.gz)")
          )
        ),
        
        tags$strong("Choose loading method:"),
        radioButtons("matrix_load_method", NULL, choices = c("Upload file" = "upload", "Load from local path" = "path"), selected = "upload", inline = TRUE),
        conditionalPanel(
          condition = "input.matrix_load_method == 'upload'",
          fileInput("scrna_matrix_file", "Upload count_matrix_sparse.mtx:", accept = c(".mtx", ".mtx.gz"))
        ),
        conditionalPanel(
          condition = "input.matrix_load_method == 'path'",
          textInput("scrna_matrix_path", "Matrix File Path:", placeholder = "/path/to/count_matrix_sparse.mtx")
        )
      ),
      
      # Column 2: Barcodes & Metadata
      card(
        card_header(class = "bg-primary text-white", "2. Genes File"),
        p("Upload the genes/features file (.tsv)."),
        
        div(style = "padding: 15px; background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 5px; margin-bottom: 15px;",
          tags$strong("Expected Format:"),
          tags$table(
            class = "table table-sm table-bordered mt-2 mb-0",
            style = "font-size: 0.85em; background: #ffffff;",
            tags$thead(
              tags$tr(style = "background: #e9ecef;",
                tags$th("Gene.symbol")
              )
            ),
            tags$tbody(
              tags$tr(tags$td("RP11-34P13.7")),
              tags$tr(tags$td("FO538757.3")),
              tags$tr(tags$td("AP006222.2")),
              tags$tr(tags$td("..."))
            )
          ),
          tags$ul(style = "margin-top: 10px; font-size: 0.9em; color: #6c757d;",
            tags$li("Single column with gene symbols"),
            tags$li("Header: 'Gene.symbol' or similar"),
            tags$li("Number of rows must match count matrix rows")
          )
        ),
        
        tags$strong("Choose loading method:"),
        radioButtons("genes_load_method", NULL, choices = c("Upload file" = "upload", "Load from local path" = "path"), selected = "upload", inline = TRUE),
        conditionalPanel(
          condition = "input.genes_load_method == 'upload'",
          fileInput("scrna_genes_file", "Upload count_matrix_genes.tsv:", accept = c(".tsv", ".txt", ".gz"))
        ),
        conditionalPanel(
          condition = "input.genes_load_method == 'path'",
          textInput("scrna_genes_path", "Genes File Path:", placeholder = "/path/to/count_matrix_genes.tsv")
        )
      )
    ),
    
    layout_column_wrap(
      width = 1/2,
      # Column 1: Barcodes
      card(
        card_header(class = "bg-primary text-white", "3. Barcodes File"),
        p("Upload the cell barcodes file (.tsv)."),
        
        div(style = "padding: 15px; background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 5px; margin-bottom: 15px;",
          tags$strong("Expected Format:"),
          tags$table(
            class = "table table-sm table-bordered mt-2 mb-0",
            style = "font-size: 0.85em; background: #ffffff;",
            tags$thead(
              tags$tr(style = "background: #e9ecef;",
                tags$th("barcode")
              )
            ),
            tags$tbody(
              tags$tr(tags$td("CID3586_AAGACCTCAGCATGAG")),
              tags$tr(tags$td("CID3586_AAGGTTCGTAGTACCT")),
              tags$tr(tags$td("CID3586_ACCAGTAGTTGTGGCC")),
              tags$tr(tags$td("..."))
            )
          ),
          tags$ul(style = "margin-top: 10px; font-size: 0.9em; color: #6c757d;",
            tags$li("Single column with cell barcodes"),
            tags$li("Header: 'barcode' or similar"),
            tags$li("Format: SampleID_BarcodeSequence, or BarcodeSequence only"),
            tags$li("Number of rows must match matrix columns"),
            tags$li("Must match barcodes in metadata file")
          )
        ),
        
        tags$strong("Choose loading method:"),
        radioButtons("barcodes_load_method", NULL, choices = c("Upload file" = "upload", "Load from local path" = "path"), selected = "upload", inline = TRUE),
        conditionalPanel(
          condition = "input.barcodes_load_method == 'upload'",
          fileInput("scrna_barcodes_file", "Upload count_matrix_barcodes.tsv:", accept = c(".tsv", ".txt", ".gz"))
        ),
        conditionalPanel(
          condition = "input.barcodes_load_method == 'path'",
          textInput("scrna_barcodes_path", "Barcodes File Path:", placeholder = "/path/to/count_matrix_barcodes.tsv")
        )
      ),
      
      # Column 2: Metadata
      card(
        card_header(class = "bg-primary text-white", "4. Cell Metadata"),
        p("Upload cell metadata with cell type annotations."),
        
        div(style = "padding: 15px; background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 5px; margin-bottom: 15px;",
          tags$strong("Expected Format:"),
          tags$div(style = "overflow-x: auto;",
            tags$table(
              class = "table table-sm table-bordered mt-2 mb-0",
              style = "font-size: 0.7em; background: #ffffff;",
              tags$thead(
                tags$tr(style = "background: #e9ecef;",
                  tags$th("barcode"),
                  tags$th("orig.ident"),
                  tags$th("nCount_RNA"),
                  tags$th("nFeature_RNA"),
                  tags$th("subtype"),
                  tags$th("celltype_major"),
                  tags$th("celltype_minor")
                )
              ),
              tags$tbody(
                tags$tr(
                  tags$td("CID3586_AACCATGAGACACC"),
                  tags$td("CID3586"),
                  tags$td("4581"),
                  tags$td("1689"),
                  tags$td("TNBC"),
                  tags$td("Endothelial"),
                  tags$td("Endothelial ACKR1")
                ),
                tags$tr(
                  tags$td("CID3586_AACCATGCATTAGC"),
                  tags$td("CID3586"),
                  tags$td("1726"),
                  tags$td("779"),
                  tags$td("TNBC"),
                  tags$td("Endothelial"),
                  tags$td("Endothelial ACKR1")
                ),
                tags$tr(
                  tags$td("CID3586_AACCATGCAGGTCACT"),
                  tags$td("CID3586"),
                  tags$td("1229"),
                  tags$td("514"),
                  tags$td("TNBC"),
                  tags$td("Endothelial"),
                  tags$td("Endothelial ACKR1")
                ),
                tags$tr(tags$td("..."), tags$td("..."), tags$td("..."), tags$td("..."), tags$td("..."), tags$td("..."), tags$td("..."))
              )
            )
          ),
          tags$div(style = "margin-top: 15px; font-size: 0.85em; color: #495057;",
            tags$strong("Key Requirements:"), tags$br(),
            tags$ul(style = "margin-top: 5px; margin-bottom: 10px; padding-left: 20px;",
              tags$li("First column: Cell barcodes (must match barcodes file exactly)"),
              tags$li("Cell type columns: Must include TWO cell type annotation columns"),
              tags$li(style = "margin-left: 20px;", "(one for crude/major cell types, one for delicate/minor cell types)"),
              tags$li("Optional columns: QC metrics (nCount_RNA, nFeature_RNA, percent.mito), sample IDs, etc.")
            ),
            tags$div(style = "margin-top: 10px;",
              tags$strong("Common delicate cell type column names:"), tags$br(),
              tags$em(style = "color: #6c757d;", "celltype_minor, celltype_subset, subcluster, annotation")
            ),
            tags$div(style = "margin-top: 8px;",
              tags$strong("Common crude cell type column names:"), tags$br(),
              tags$em(style = "color: #6c757d;", "celltype_major, cell_type, CellType, cluster")
            )
          )
        ),
        tags$div(style = "font-size: 0.9em; color: #6c757d; margin-bottom: 15px; font-style: italic;",
          "Required: Must contain cell type annotations. Supported formats: .csv, .tsv, .txt"
        ),
        
        tags$strong("Choose loading method:"),
        radioButtons("metadata_load_method", NULL, choices = c("Upload file" = "upload", "Load from local path" = "path"), selected = "upload", inline = TRUE),
        conditionalPanel(
          condition = "input.metadata_load_method == 'upload'",
          fileInput("scrna_metadata_file", "Upload metadata.csv:", accept = c(".csv", ".tsv"))
        ),
        conditionalPanel(
          condition = "input.metadata_load_method == 'path'",
          textInput("scrna_metadata_path", "Metadata File Path:", placeholder = "/path/to/metadata.csv")
        )
      )
    ),
    
    card(
      actionButton("load_scrna_btn", "Load Raw Data", class = "btn btn-primary btn-lg"),
      hr(),
      uiOutput("scrna_files_status"),
      uiOutput("scrna_upload_status")
    )
  ),

  # --- Step 4: scRNA-seq Processing ---
  nav_panel(
    title = "Step 4: Processing", 
    value = "scrna_processing",
    icon = bs_icon("gear"),
    card(
      card_header(class = "bg-primary text-white", "scRNA-seq Data Processing & Validation"),
      
      layout_column_wrap(
        width = 1/3,
        card(
          card_header("Sample Selection (Optional)"),
          p("Filter cells by sample subtype or patient resource. Leave blank to keep all samples."),
          selectInput("subtype_column", "Subtype Column Name:", choices = NULL),
          tags$div(style = "font-size: 0.9em; color: #6c757d; font-style: italic; margin-bottom: 10px;",
            "Select the metadata column containing sample subtypes and choose which subtypes to keep."
          ),
          uiOutput("subtype_values_ui")
        ),
        
        card(
          card_header("Gene Filtering"),
          numericInput("min_cells_per_gene", "Minimum cells expressing each gene:", value = 10, min = 1, max = 100),
          tags$div(style = "font-size: 0.9em; color: #6c757d; font-style: italic;",
            "Genes expressed in fewer cells will be removed."
          )
        ),
        
        card(
          card_header("Cell Type Annotation Columns"),
          p("Specify the two cell type annotation columns (major and minor)."),
          tags$strong("Major Cell Type Column:"),
          selectInput("celltype_major_column", NULL, choices = NULL),
          tags$strong("Minor Cell Type Column:"),
          selectInput("celltype_minor_column", NULL, choices = NULL),
          tags$div(style = "font-size: 0.85em; color: #6c757d; font-style: italic; margin-top: 10px;",
            "Major = broader categories (e.g., T cells, B cells). Minor = specific subtypes (e.g., CD4+ T cells, CD8+ T cells)."
          )
        )
      ),
      
      card(
        card_header("Run Processing"),
        actionButton("process_scrna_btn", "Run Preprocessing & Validation", class = "btn btn-primary w-100 btn-lg", icon = bs_icon("gear"))
      ),
      
      hr(),
      tags$h5("Processing Results", style = "color: #002F6C;"),
      verbatimTextOutput("scrna_processing_summary")
    )
  ),

  # --- Step 5: Query Method ---
  nav_panel(
    title = "Step 5: Query Method", 
    value = "query_method",
    icon = bs_icon("search"),
    card(
      card_header(class = "bg-primary text-white", "Cell-Type Enrichment Analysis"),
      div(
        class = "alert alert-info",
        p(bs_icon("info-circle"), " The Query Method computes signature scores for region-specific proteins across cell types, identifying enriched populations using z-score statistics.")
      ),
      actionButton("run_query_method_btn", "Run Query Method Analysis", class = "btn btn-primary btn-lg"),
      hr(),
      uiOutput("query_method_status")
    ),
    navset_card_pill(
      nav_panel(
        "Positive Enrichment",
        downloadButton("download_query_plot", "Download Plot", class = "btn-sm mb-3"),
        plotOutput("query_bubble_plot_positive", height = "600px")
      ),
      nav_panel(
        "Negative Enrichment",
        downloadButton("download_query_plot_negative", "Download Plot", class = "btn-sm mb-3"),
        plotOutput("query_bubble_plot_negative", height = "600px")
      ),
      nav_panel(
        "Overlap Analysis",
        downloadButton("download_overlap_genes", "Download Overlap Genes", class = "btn-sm mb-3"),
        DTOutput("overlap_summary_table"),
        hr(),
        DTOutput("overlap_genes_table")
      ),
      nav_panel(
        "Results Table",
        downloadButton("download_query_results", "Download Results", class = "btn-sm mb-3"),
        DTOutput("query_results_table")
      )
    )
  ),

  # --- Results Tab ---
  nav_panel(
    title = "Results & Download", 
    value = "results",
    icon = bs_icon("download"),
    card(
      card_header(class = "bg-primary text-white", "Analysis Summary & Downloads")
    ),
    card(
      card_header("Analysis Summary"),
      verbatimTextOutput("analysis_summary")
    ),
    layout_column_wrap(
      width = 1/3,
      card(
        card_header("Download All Results"),
        p("Download a ZIP file containing all analysis results, tables, and plots."),
        downloadButton("download_all_results", "Download ZIP", class = "btn btn-primary w-100")
      ),
      card(
        card_header("Generate Report"),
        p("Generate a comprehensive PDF report of the analysis."),
        downloadButton("download_report", "Generate PDF Report", class = "btn btn-primary w-100")
      ),
      card(
        card_header("Reset Analysis"),
        p(class = "text-danger", "Warning: This will clear all data and results."),
        actionButton("reset_analysis", "Reset Analysis", class = "btn btn-primary w-100")
      )
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  
  # File upload size limit is set at the top of the file (10GB)
  # No need to set it again here
  
  # Reactive values to store data
  rv <- reactiveValues(
    spatial_expr = NULL,
    metadata = NULL,
    scrna_obj = NULL,
    scrna_obj_filtered = NULL,
    scrna_raw_matrix = NULL,
    scrna_raw_metadata = NULL,
    filtered_genes = NULL,
    region_specific_genes = NULL,
    region_gene_list = NULL,
    region_gene_list_overlap = NULL,
    query_results = NULL,
    pca_result = NULL,
    pca_data = NULL,
    var_explained = NULL,
    overlap_matrix = NULL,
    region_gene_list_for_venn = NULL,
    # scRNA-seq file status tracking
    scrna_matrix_loaded = FALSE,
    scrna_genes_loaded = FALSE,
    scrna_barcodes_loaded = FALSE,
    scrna_metadata_loaded = FALSE,
    scrna_matrix_data = NULL,
    scrna_genes_data = NULL,
    scrna_barcodes_data = NULL,
    scrna_metadata_data = NULL,
    proteomics_validated = FALSE,
    scrna_loaded = FALSE,
    scrna_processed = FALSE,
    region_analysis_complete = FALSE,
    query_method_complete = FALSE
  )
  
  # Navigate to upload page when "Get Started" is clicked
  observeEvent(input$start_btn, {
    nav_select("navbar", "proteomics_upload")
  })
  
  # Load Spatial Expression Data
  observeEvent(input$spatial_expr_file, {
    req(input$spatial_expr_file)
    
    tryCatch({
      file_path <- input$spatial_expr_file$datapath
      file_ext <- tools::file_ext(input$spatial_expr_file$name)
      
      if (file_ext %in% c("xlsx", "xls")) {
        rv$spatial_expr <- read_excel(file_path)
      } else if (file_ext == "csv") {
        rv$spatial_expr <- read.csv(file_path, check.names = FALSE)
      } else if (file_ext %in% c("txt", "tsv")) {
        rv$spatial_expr <- read.delim(file_path, check.names = FALSE)
      }
      
      # Basic validation
      n_genes <- nrow(rv$spatial_expr)
      n_samples <- ncol(rv$spatial_expr) - 1
      gene_col_name <- colnames(rv$spatial_expr)[1]
      sample_ids <- colnames(rv$spatial_expr)[-1]
      
      output$spatial_expr_preview <- renderUI({
        div(
          style = "background-color: #dff0d8; padding: 10px; border-radius: 5px;",
          p(icon("check-circle"), tags$b("✓ Expression data loaded successfully!"), style = "color: green; font-weight: bold;"),
          p(paste("📊 Dimensions:", n_genes, "genes ×", n_samples, "samples")),
          p(paste("🧬 Gene column:", gene_col_name)),
          p(paste("📝 Sample IDs:", paste(head(sample_ids, 5), collapse = ", "), 
                  ifelse(length(sample_ids) > 5, "...", ""))),
          div(
            style = "margin-top: 10px; padding: 8px; background-color: #d9edf7; border-radius: 3px;",
            p(style = "margin: 0; font-size: 12px;",
              icon("info-circle"), " Preview: ", nrow(rv$spatial_expr), " genes, ",
              "First gene: ", rv$spatial_expr[1,1], ", ",
              round(max(rv$spatial_expr[,-1], na.rm=TRUE), 2))
          )
        )
      })
    }, error = function(e) {
      output$spatial_expr_preview <- renderUI({
        div(
          style = "background-color: #f2dede; padding: 10px; border-radius: 5px;",
          p(icon("exclamation-circle"), tags$b("Error loading file:"), style = "color: red; font-weight: bold;"),
          p(e$message, style = "color: red;")
        )
      })
    })
  })
  
  # Load Metadata
  observeEvent(input$metadata_file, {
    req(input$metadata_file)
    
    tryCatch({
      file_path <- input$metadata_file$datapath
      file_ext <- tools::file_ext(input$metadata_file$name)
      
      if (file_ext %in% c("xlsx", "xls")) {
        rv$metadata <- read_excel(file_path)
      } else if (file_ext == "csv") {
        rv$metadata <- read.csv(file_path, check.names = FALSE)
      } else if (file_ext %in% c("txt", "tsv")) {
        rv$metadata <- read.delim(file_path, check.names = FALSE)
      }
      
      # Update region column choices
      updateSelectInput(session, "group_column", choices = colnames(rv$metadata))
      
      # Get sample info
      n_samples <- nrow(rv$metadata)
      sample_col <- colnames(rv$metadata)[1]
      
      # Try to find region/group column
      region_col <- grep("group|region|type|condition", colnames(rv$metadata), 
                         ignore.case = TRUE, value = TRUE)[1]
      if (is.na(region_col)) region_col <- colnames(rv$metadata)[2]
      
      # Set default selection
      updateSelectInput(session, "group_column", selected = region_col)
      
      unique_groups <- unique(rv$metadata[[region_col]])
      
      output$metadata_preview <- renderUI({
        msg_list <- list(
          div(
            style = "background-color: #dff0d8; padding: 10px; border-radius: 5px;",
            p(icon("check-circle"), tags$b("✓ Metadata loaded successfully!"), style = "color: green; font-weight: bold;"),
            p(paste("📋 Samples:", n_samples)),
            p(paste("📊 Columns:", paste(colnames(rv$metadata), collapse = ", "))),
            p(paste("🏷️ Detected region column:", region_col)),
            p(paste("🔢 Unique groups:", length(unique_groups), "-", 
                    paste(head(unique_groups, 3), collapse = ", "),
                    ifelse(length(unique_groups) > 3, "...", "")))
          )
        )
        
        # Check if expression data is loaded and samples match
        if (!is.null(rv$spatial_expr)) {
          expr_samples <- colnames(rv$spatial_expr)[-1]
          meta_samples <- rv$metadata[[sample_col]]
          matching_samples <- intersect(expr_samples, meta_samples)
          
          if (length(matching_samples) > 0) {
            msg_list[[2]] <- div(
              style = "margin-top: 10px; padding: 8px; background-color: #d9edf7; border-radius: 3px;",
              p(style = "margin: 0; font-size: 12px; color: #31708f;",
                icon("check"), " Sample ID Match: ", length(matching_samples), "/", 
                length(expr_samples), " samples match expression data")
            )
          } else {
            msg_list[[2]] <- div(
              style = "margin-top: 10px; padding: 8px; background-color: #fcf8e3; border-radius: 3px;",
              p(style = "margin: 0; font-size: 12px; color: #8a6d3b;",
                icon("warning"), " Warning: Sample IDs don't match expression data column names!")
            )
          }
        }
        
        tagList(msg_list)
      })
    }, error = function(e) {
      output$metadata_preview <- renderUI({
        div(
          style = "background-color: #f2dede; padding: 10px; border-radius: 5px;",
          p(icon("exclamation-circle"), tags$b("Error loading file:"), style = "color: red; font-weight: bold;"),
          p(e$message, style = "color: red;")
        )
      })
    })
  })
  
  # Track scRNA-seq file uploads and paths
  
  # Matrix file upload
  observeEvent(input$scrna_matrix_file, {
    req(input$scrna_matrix_file)
    req(input$matrix_load_method == "upload")
    
    tryCatch({
      rv$scrna_matrix_data <- input$scrna_matrix_file$datapath
      rv$scrna_matrix_loaded <- TRUE
      showNotification("✓ Matrix file uploaded", type = "message", duration = 3)
      update_scrna_files_status()
    }, error = function(e) {
      rv$scrna_matrix_loaded <- FALSE
      showNotification(paste("Error uploading matrix file:", e$message), type = "error", duration = 5)
    })
  })
  
  # Matrix file path
  observeEvent(input$scrna_matrix_path, {
    req(input$scrna_matrix_path)
    req(input$matrix_load_method == "path")
    
    file_path <- trimws(input$scrna_matrix_path)
    if (file_path != "" && file.exists(file_path)) {
      rv$scrna_matrix_data <- file_path
      rv$scrna_matrix_loaded <- TRUE
      showNotification("✓ Matrix file path validated", type = "message", duration = 3)
      update_scrna_files_status()
    } else if (file_path != "") {
      rv$scrna_matrix_loaded <- FALSE
      showNotification("File not found. Check the path.", type = "warning", duration = 5)
      update_scrna_files_status()
    }
  })
  
  # Genes file upload
  observeEvent(input$scrna_genes_file, {
    req(input$scrna_genes_file)
    req(input$genes_load_method == "upload")
    
    tryCatch({
      rv$scrna_genes_data <- input$scrna_genes_file$datapath
      rv$scrna_genes_loaded <- TRUE
      showNotification("✓ Genes file uploaded", type = "message", duration = 3)
      update_scrna_files_status()
    }, error = function(e) {
      rv$scrna_genes_loaded <- FALSE
      showNotification(paste("Error uploading genes file:", e$message), type = "error", duration = 5)
    })
  })
  
  # Genes file path
  observeEvent(input$scrna_genes_path, {
    req(input$scrna_genes_path)
    req(input$genes_load_method == "path")
    
    file_path <- trimws(input$scrna_genes_path)
    if (file_path != "" && file.exists(file_path)) {
      rv$scrna_genes_data <- file_path
      rv$scrna_genes_loaded <- TRUE
      showNotification("✓ Genes file path validated", type = "message", duration = 3)
      update_scrna_files_status()
    } else if (file_path != "") {
      rv$scrna_genes_loaded <- FALSE
      showNotification("File not found. Check the path.", type = "warning", duration = 5)
      update_scrna_files_status()
    }
  })
  
  # Barcodes file upload
  observeEvent(input$scrna_barcodes_file, {
    req(input$scrna_barcodes_file)
    req(input$barcodes_load_method == "upload")
    
    tryCatch({
      rv$scrna_barcodes_data <- input$scrna_barcodes_file$datapath
      rv$scrna_barcodes_loaded <- TRUE
      showNotification("✓ Barcodes file uploaded", type = "message", duration = 3)
      update_scrna_files_status()
    }, error = function(e) {
      rv$scrna_barcodes_loaded <- FALSE
      showNotification(paste("Error uploading barcodes file:", e$message), type = "error", duration = 5)
    })
  })
  
  # Barcodes file path
  observeEvent(input$scrna_barcodes_path, {
    req(input$scrna_barcodes_path)
    req(input$barcodes_load_method == "path")
    
    file_path <- trimws(input$scrna_barcodes_path)
    if (file_path != "" && file.exists(file_path)) {
      rv$scrna_barcodes_data <- file_path
      rv$scrna_barcodes_loaded <- TRUE
      showNotification("✓ Barcodes file path validated", type = "message", duration = 3)
      update_scrna_files_status()
    } else if (file_path != "") {
      rv$scrna_barcodes_loaded <- FALSE
      showNotification("File not found. Check the path.", type = "warning", duration = 5)
      update_scrna_files_status()
    }
  })
  
  # Metadata file upload
  observeEvent(input$scrna_metadata_file, {
    req(input$scrna_metadata_file)
    req(input$metadata_load_method == "upload")
    
    tryCatch({
      rv$scrna_metadata_data <- input$scrna_metadata_file$datapath
      rv$scrna_metadata_loaded <- TRUE
      showNotification("✓ Metadata file uploaded", type = "message", duration = 3)
      update_scrna_files_status()
    }, error = function(e) {
      rv$scrna_metadata_loaded <- FALSE
      showNotification(paste("Error uploading metadata file:", e$message), type = "error", duration = 5)
    })
  })
  
  # Metadata file path
  observeEvent(input$scrna_metadata_path, {
    req(input$scrna_metadata_path)
    req(input$metadata_load_method == "path")
    
    file_path <- trimws(input$scrna_metadata_path)
    if (file_path != "" && file.exists(file_path)) {
      rv$scrna_metadata_data <- file_path
      rv$scrna_metadata_loaded <- TRUE
      showNotification("✓ Metadata file path validated", type = "message", duration = 3)
      update_scrna_files_status()
    } else if (file_path != "") {
      rv$scrna_metadata_loaded <- FALSE
      showNotification("File not found. Check the path.", type = "warning", duration = 5)
      update_scrna_files_status()
    }
  })
  
  # Update file upload status display
  update_scrna_files_status <- function() {
    output$scrna_files_status <- renderUI({
      div(
        style = "background-color: #f9f9f9; padding: 15px; border-radius: 5px; margin: 10px 0;",
        h4("File Upload Status:"),
        tags$ul(
          tags$li(
            if (rv$scrna_matrix_loaded) {
              tags$span(icon("check-circle"), " Matrix file", style = "color: green; font-weight: bold;")
            } else {
              tags$span(icon("circle"), " Matrix file", style = "color: gray;")
            }
          ),
          tags$li(
            if (rv$scrna_genes_loaded) {
              tags$span(icon("check-circle"), " Genes file", style = "color: green; font-weight: bold;")
            } else {
              tags$span(icon("circle"), " Genes file", style = "color: gray;")
            }
          ),
          tags$li(
            if (rv$scrna_barcodes_loaded) {
              tags$span(icon("check-circle"), " Barcodes file", style = "color: green; font-weight: bold;")
            } else {
              tags$span(icon("circle"), " Barcodes file", style = "color: gray;")
            }
          ),
          tags$li(
            if (rv$scrna_metadata_loaded) {
              tags$span(icon("check-circle"), " Metadata file", style = "color: green; font-weight: bold;")
            } else {
              tags$span(icon("circle"), " Metadata file", style = "color: gray;")
            }
          )
        ),
        if (rv$scrna_matrix_loaded && rv$scrna_genes_loaded && 
            rv$scrna_barcodes_loaded && rv$scrna_metadata_loaded) {
          div(
            style = "margin-top: 10px; padding: 10px; background-color: #dff0d8; border-radius: 3px;",
            p(icon("check"), tags$b(" All files uploaded! Ready to create Seurat object."),
              style = "color: green; margin: 0;")
          )
        } else {
          div(
            style = "margin-top: 10px; padding: 10px; background-color: #fcf8e3; border-radius: 3px;",
            p(icon("info-circle"), " Please upload all 4 files before proceeding.",
              style = "color: #8a6d3b; margin: 0;")
          )
        }
      )
    })
  }
  
  # Initialize file status display
  output$scrna_files_status <- renderUI({
    div(
      style = "background-color: #f9f9f9; padding: 15px; border-radius: 5px; margin: 10px 0;",
      h4("File Upload Status:"),
      tags$ul(
        tags$li(tags$span(icon("circle"), " Matrix file", style = "color: gray;")),
        tags$li(tags$span(icon("circle"), " Genes file", style = "color: gray;")),
        tags$li(tags$span(icon("circle"), " Barcodes file", style = "color: gray;")),
        tags$li(tags$span(icon("circle"), " Metadata file", style = "color: gray;"))
      ),
      div(
        style = "margin-top: 10px; padding: 10px; background-color: #fcf8e3; border-radius: 3px;",
        p(icon("info-circle"), " Please upload all 4 files before proceeding.",
          style = "color: #8a6d3b; margin: 0;")
      )
    )
  })
  
  # Validate Proteomics Data (Step 1)
  observeEvent(input$validate_proteomics_btn, {
    validation_messages <- c()
    errors <- c()
    warnings <- c()
    
    # ======================
    # 1. Check if files are loaded
    # ======================
    if (is.null(rv$spatial_expr)) {
      errors <- c(errors, "✗ Spatial expression data not loaded")
    } else {
      validation_messages <- c(validation_messages, "✓ Spatial expression data loaded")
      
      # Check expression data structure
      if (ncol(rv$spatial_expr) < 2) {
        errors <- c(errors, "  ✗ Expression data must have at least 2 columns (genes + samples)")
      } else {
        validation_messages <- c(validation_messages, 
                                 paste("  • Dimensions:", nrow(rv$spatial_expr), "genes ×", 
                                       ncol(rv$spatial_expr) - 1, "samples"))
      }
      
      # Check for numeric values
      numeric_cols <- sapply(rv$spatial_expr[, -1], is.numeric)
      if (!all(numeric_cols)) {
        errors <- c(errors, paste("  ✗ Non-numeric values found in", sum(!numeric_cols), "sample column(s)"))
      }
      
      # Check for missing values
      na_count <- sum(is.na(rv$spatial_expr[, -1]))
      if (na_count > 0) {
        warnings <- c(warnings, paste("  ⚠ Warning:", na_count, "missing values in expression data"))
      }
      
      # Check for duplicate gene names
      gene_col <- rv$spatial_expr[[1]]
      dup_genes <- sum(duplicated(gene_col))
      if (dup_genes > 0) {
        errors <- c(errors, paste("  ✗", dup_genes, "duplicate gene names found"))
      }
    }
    
    if (is.null(rv$metadata)) {
      errors <- c(errors, "✗ Metadata not loaded")
    } else {
      validation_messages <- c(validation_messages, "✓ Metadata loaded")
      
      # Check metadata structure
      if (ncol(rv$metadata) < 2) {
        errors <- c(errors, "  ✗ Metadata must have at least 2 columns (Sample_ID + Group)")
      } else {
        validation_messages <- c(validation_messages, 
                                 paste("  • Samples:", nrow(rv$metadata),
                                       "| Columns:", ncol(rv$metadata)))
      }
      
      # Check for missing values in metadata
      na_count_meta <- sum(is.na(rv$metadata))
      if (na_count_meta > 0) {
        warnings <- c(warnings, paste("  ⚠ Warning:", na_count_meta, "missing values in metadata"))
      }
      
      # Check for group column
      if (ncol(rv$metadata) >= 2) {
        region_col <- grep("group|region|type|condition", colnames(rv$metadata), 
                           ignore.case = TRUE, value = TRUE)[1]
        if (is.na(region_col)) region_col <- colnames(rv$metadata)[2]
        
        unique_groups <- unique(rv$metadata[[region_col]])
        unique_groups <- unique_groups[!is.na(unique_groups)]
        
        if (length(unique_groups) < 2) {
          errors <- c(errors, "  ✗ At least 2 different groups/regions required for comparison")
        } else {
          validation_messages <- c(validation_messages,
                                   paste("  • Detected region column:", region_col),
                                   paste("  • Unique groups:", length(unique_groups)))
        }
      }
    }
    
    # ======================
    # 2. Cross-validation between files
    # ======================
    if (!is.null(rv$spatial_expr) && !is.null(rv$metadata)) {
      validation_messages <- c(validation_messages, "", "═══ Cross-Validation ═══")
      
      # Get sample IDs
      expr_samples <- colnames(rv$spatial_expr)[-1]
      meta_samples <- as.character(rv$metadata[[1]])
      
      # Check sample ID matching
      matching_samples <- intersect(expr_samples, meta_samples)
      missing_in_meta <- setdiff(expr_samples, meta_samples)
      missing_in_expr <- setdiff(meta_samples, expr_samples)
      
      if (length(matching_samples) == 0) {
        errors <- c(errors, "✗ CRITICAL: No matching sample IDs between expression and metadata!")
        errors <- c(errors, "  • Check that sample names match exactly (case-sensitive)")
      } else if (length(matching_samples) == length(expr_samples) && 
                 length(matching_samples) == length(meta_samples)) {
        validation_messages <- c(validation_messages,
                                 paste("✓ Perfect match:", length(matching_samples), "samples"))
      } else {
        validation_messages <- c(validation_messages,
                                 paste("✓ Partial match:", length(matching_samples), "samples"))
        
        if (length(missing_in_meta) > 0) {
          warnings <- c(warnings, 
                        paste("  ⚠", length(missing_in_meta), 
                              "sample(s) in expression not found in metadata:",
                              paste(head(missing_in_meta, 3), collapse = ", "),
                              ifelse(length(missing_in_meta) > 3, "...", "")))
        }
        
        if (length(missing_in_expr) > 0) {
          warnings <- c(warnings,
                        paste("  ⚠", length(missing_in_expr),
                              "sample(s) in metadata not found in expression:",
                              paste(head(missing_in_expr, 3), collapse = ", "),
                              ifelse(length(missing_in_expr) > 3, "...", "")))
        }
      }
    }
    
    # ======================
    # 3. Final validation result
    # ======================
    rv$proteomics_validated <- (length(errors) == 0 && 
                                  !is.null(rv$spatial_expr) && 
                                  !is.null(rv$metadata))
    
    # Combine all messages
    all_messages <- c(validation_messages)
    
    if (length(warnings) > 0) {
      all_messages <- c(all_messages, "", "═══ Warnings ═══", warnings)
    }
    
    if (length(errors) > 0) {
      all_messages <- c(all_messages, "", "═══ Errors ═══", errors)
    }
    
    if (rv$proteomics_validated) {
      all_messages <- c(all_messages, "", 
                        "════════════════════════════════════",
                        "✓✓✓ All data validated successfully!",
                        "You can proceed to Step 2: Region-Specific Proteins",
                        "════════════════════════════════════")
      showNotification("✓ Data validation successful!", type = "message", duration = 5)
    } else if (length(errors) > 0) {
      all_messages <- c(all_messages, "", 
                        "════════════════════════════════════",
                        "✗ Validation failed! Please fix the errors above.",
                        "════════════════════════════════════")
      showNotification("✗ Validation failed. Check errors above.", type = "error", duration = 10)
    } else {
      all_messages <- c(all_messages, "", 
                        "════════════════════════════════════",
                        "⚠ Please upload all required files before proceeding.",
                        "════════════════════════════════════")
      showNotification("⚠ Please upload all required files", type = "warning", duration = 5)
    }
    
    output$proteomics_validation_summary <- renderText({
      paste(all_messages, collapse = "\n")
    })
  })
  
  # Load scRNA-seq Data (Step 3)
  # Load scRNA-seq Data (Step 3) - Create Seurat object from 10X files
  observeEvent(input$load_scrna_btn, {
    
    # Check if all 4 files are uploaded
    if (!rv$scrna_matrix_loaded || !rv$scrna_genes_loaded || 
        !rv$scrna_barcodes_loaded || !rv$scrna_metadata_loaded) {
      showNotification("⚠ Please upload all 4 required files before proceeding.", 
                       type = "warning", duration = 5)
      return()
    }
    
    withProgress(message = 'Creating Seurat object from 10X data...', value = 0, {
      
      tryCatch({
        incProgress(0.1, detail = "Reading matrix file...")
        
        # Read sparse matrix
        count_matrix <- readMM(rv$scrna_matrix_data)
        
        incProgress(0.3, detail = "Reading genes file...")
        
        # Read genes file - use scan() with skip=1 for TSV with header
        gene_names <- scan(rv$scrna_genes_data, character(), skip = 1)
        gene_names <- trimws(gene_names)
        
        incProgress(0.4, detail = "Reading barcodes file...")
        
        # Read barcodes file - use scan() with skip=1 for TSV with header
        barcode_names <- scan(rv$scrna_barcodes_data, character(), skip = 1)
        barcode_names <- trimws(barcode_names)
        
        incProgress(0.5, detail = "Reading metadata file...")
        
        # Read metadata file
        file_ext <- tools::file_ext(input$scrna_metadata_file$name)
        if (file_ext == "csv") {
          metadata <- read.csv(rv$scrna_metadata_data, stringsAsFactors = FALSE)
        } else {
          metadata <- read.delim(rv$scrna_metadata_data, stringsAsFactors = FALSE)
        }
        
        incProgress(0.6, detail = "Validating data dimensions...")
        
        # Validate dimensions
        if (nrow(count_matrix) != length(gene_names)) {
          stop(paste("Dimension mismatch: Matrix has", nrow(count_matrix), 
                     "rows but genes file has", length(gene_names), "entries"))
        }
        
        if (ncol(count_matrix) != length(barcode_names)) {
          stop(paste("Dimension mismatch: Matrix has", ncol(count_matrix), 
                     "columns but barcodes file has", length(barcode_names), "entries"))
        }
        
        # Check if metadata barcodes match matrix barcodes
        meta_barcode_col <- colnames(metadata)[1]
        meta_barcodes <- metadata[[meta_barcode_col]]
        
        matching_barcodes <- intersect(barcode_names, meta_barcodes)
        
        if (length(matching_barcodes) == 0) {
          stop("No matching barcodes between matrix and metadata. Please check that barcode formats match.")
        }
        
        if (length(matching_barcodes) < length(barcode_names)) {
          warning_msg <- paste(length(barcode_names) - length(matching_barcodes), 
                               "barcodes in matrix not found in metadata. These cells will be excluded.")
        }
        
        incProgress(0.7, detail = "Setting matrix dimensions...")
        
        # Set row and column names
        rownames(count_matrix) <- gene_names
        colnames(count_matrix) <- barcode_names
        
        # Ensure gene names are unique to avoid Seurat LogMap errors
        if (any(duplicated(rownames(count_matrix)))) {
          rownames(count_matrix) <- make.unique(rownames(count_matrix))
        }
        
        # Reorder count_matrix columns to match metadata order (like user's match approach)
        # User's code: data <- data[,match(metadata[[1]], colnames(data))]
        count_matrix <- count_matrix[, match(metadata[[1]], colnames(count_matrix))]
        
        incProgress(0.8, detail = "Storing raw data...")
        
        # Store raw data without creating Seurat object yet
        # This allows filtering by subtype in Step 4 before object creation
        # DO NOT set rownames on metadata here - will be set in Step 4 after filtering
        rv$scrna_raw_matrix <- count_matrix
        rv$scrna_raw_metadata <- metadata
        
        rv$scrna_loaded <- TRUE
        
        incProgress(1, detail = "Complete!")
        
        # Display success message
        n_cells <- ncol(rv$scrna_raw_matrix)
        n_genes <- nrow(rv$scrna_raw_matrix)
        n_meta_cols <- ncol(rv$scrna_raw_metadata)
        
        output$scrna_upload_status <- renderUI({
          div(
            style = "background-color: #dff0d8; padding: 15px; border-radius: 5px; margin: 10px 0;",
            h4(icon("check-circle"), "scRNA-seq Data Loaded Successfully!", style = "color: green;"),
            p(paste("📊 Total Cells:", n_cells)),
            p(paste("🧬 Total Genes:", n_genes)),
            p(paste("📋 Metadata columns:", n_meta_cols)),
            p(paste("🔖 Available metadata:", paste(head(colnames(rv$scrna_raw_metadata), 5), collapse = ", "),
                    ifelse(n_meta_cols > 5, "...", ""))),
            hr(),
            div(
              style = "background-color: #d9edf7; padding: 10px; border-radius: 3px;",
              p(icon("arrow-right"), tags$b(" Next Step:"), " Proceed to Step 4 to filter by subtype (optional) and create Seurat object.",
                style = "margin: 0;")
            ),
            div(
              style = "background-color: #fcf8e3; padding: 10px; border-radius: 3px; margin-top: 10px;",
              p(icon("lightbulb"), tags$b(" Performance Tip:"), " Filtering by subtype before creating Seurat object can significantly reduce computation time!",
                style = "margin: 0;")
            )
          )
        })
        
        showNotification("✓ scRNA-seq data loaded successfully!", type = "message", duration = 5)
        
        # Update Step 4 dropdown menus with metadata columns
        metadata_cols <- colnames(rv$scrna_raw_metadata)
        
        # Try to auto-detect subtype column
        subtype_col <- grep("subtype|type|condition|group", metadata_cols, 
                            ignore.case = TRUE, value = TRUE)[1]
        
        updateSelectInput(session, "subtype_column",
                          choices = c("None (keep all)" = "none", metadata_cols),
                          selected = if(!is.na(subtype_col)) subtype_col else "none")
        
        # Try to auto-detect cell type columns
        major_col <- grep("celltype_major|cell_type|CellType|cluster", metadata_cols, 
                          ignore.case = TRUE, value = TRUE)[1]
        minor_col <- grep("celltype_minor|celltype_subset|subcluster|annotation", metadata_cols,
                          ignore.case = TRUE, value = TRUE)[1]
        
        updateSelectInput(session, "celltype_major_column",
                          choices = metadata_cols,
                          selected = if(!is.na(major_col)) major_col else metadata_cols[1])
        
        updateSelectInput(session, "celltype_minor_column",
                          choices = metadata_cols,
                          selected = if(!is.na(minor_col)) minor_col else metadata_cols[1])
        
      }, error = function(e) {
        output$scrna_upload_status <- renderUI({
          div(
            style = "background-color: #f2dede; padding: 15px; border-radius: 5px;",
            h4(icon("exclamation-circle"), "Error Loading scRNA-seq Data:", style = "color: red;"),
            p(e$message, style = "color: red;"),
            hr(),
            div(
              style = "background-color: #fcf8e3; padding: 10px; border-radius: 3px;",
              p(tags$b("Troubleshooting Tips:"), style = "margin: 0 0 5px 0;"),
              tags$ul(
                style = "margin: 5px 0; padding-left: 20px;",
                tags$li("Ensure matrix dimensions match genes and barcodes files"),
                tags$li("Check that barcode format in metadata matches barcodes file"),
                tags$li("Verify all files are from the same 10X dataset"),
                tags$li("Check for any non-ASCII characters in files")
              )
            )
          )
        })
        showNotification(paste("✗ Error:", e$message), type = "error", duration = 10)
      })
    })
  })
  
  # Update subtype values when subtype column is selected
  observeEvent(input$subtype_column, {
    req(rv$scrna_loaded)
    req(input$subtype_column)
    
    if (input$subtype_column != "none") {
      subtype_values <- unique(rv$scrna_raw_metadata[[input$subtype_column]])
      subtype_values <- subtype_values[!is.na(subtype_values)]
      
      output$subtype_values_ui <- renderUI({
        div(
          selectInput("subtype_values",
                      "Select Subtypes to Keep:",
                      choices = subtype_values,
                      selected = subtype_values,
                      multiple = TRUE),
          p(class = "example-text", style = "color: #3c8dbc;",
            paste("💡 Filtering will reduce", ncol(rv$scrna_raw_matrix), "cells to selected subtypes before creating Seurat object."))
        )
      })
    } else {
      output$subtype_values_ui <- renderUI({
        p(class = "example-text", "All samples will be kept.")
      })
    }
  })
  
  # Process & Filter scRNA-seq Data (Step 4)
  observeEvent(input$process_scrna_btn, {
    req(rv$scrna_loaded)
    req(rv$scrna_raw_matrix)
    req(rv$scrna_raw_metadata)
    
    withProgress(message = 'Preprocessing scRNA-seq data...', value = 0, {
      
      tryCatch({
        incProgress(0.1, detail = "Starting preprocessing...")
        
        processing_log <- c()
        
        # Start with raw data
        count_matrix <- rv$scrna_raw_matrix
        metadata <- rv$scrna_raw_metadata
        
        processing_log <- c(processing_log, paste("Initial data:", ncol(count_matrix), "cells,", nrow(count_matrix), "genes"))
        
        # Filter by subtype if selected
        incProgress(0.2, detail = "Filtering by subtype...")
        if (!is.null(input$subtype_column) && input$subtype_column != "none" && !is.null(input$subtype_values)) {
          subtype_col <- input$subtype_column
          selected_subtypes <- input$subtype_values
          
          # Get indices of cells matching selected subtypes (like user's grep approach)
          keep_inds <- which(metadata[[subtype_col]] %in% selected_subtypes)
          
          if (length(keep_inds) == 0) {
            stop("No cells found with selected subtypes")
          }
          
          # Filter both metadata and count_matrix by the same indices
          metadata <- metadata[keep_inds, ]
          count_matrix <- count_matrix[, keep_inds]
          
          processing_log <- c(processing_log, 
                              paste("Filtered by", subtype_col, "to keep:", paste(selected_subtypes, collapse = ", ")),
                              paste("After subtype filtering:", ncol(count_matrix), "cells remain"))
        } else {
          processing_log <- c(processing_log, "No subtype filtering applied - keeping all cells")
        }
        
        # Set rownames AFTER filtering (like user's code: rownames(meta.combine) <- meta.combine[[1]])
        barcode_ids <- metadata[[1]]
        if (any(duplicated(barcode_ids))) {
          dup_count <- sum(duplicated(barcode_ids))
          processing_log <- c(processing_log, paste("Warning:", dup_count, "duplicate barcodes in metadata. Making unique..."))
          barcode_ids <- make.unique(as.character(barcode_ids))
        }
        
        if (any(duplicated(colnames(count_matrix)))) {
          dup_count <- sum(duplicated(colnames(count_matrix)))
          processing_log <- c(processing_log, paste("Warning:", dup_count, "duplicate barcodes in count matrix. Making unique..."))
          colnames(count_matrix) <- make.unique(colnames(count_matrix))
        }
        
        if (any(is.na(barcode_ids)) || any(barcode_ids == "")) {
          stop("ERROR: Empty or NA barcodes found in metadata after filtering")
        }
        
        if (any(is.na(rownames(count_matrix))) || any(rownames(count_matrix) == "")) {
          stop("ERROR: Empty or NA gene names found in count matrix")
        }
        
        if (any(duplicated(rownames(count_matrix)))) {
          dup_count <- sum(duplicated(rownames(count_matrix)))
          processing_log <- c(processing_log, paste("Warning:", dup_count, "duplicate gene names in count matrix. Making unique..."))
          rownames(count_matrix) <- make.unique(rownames(count_matrix))
        }
        
        rownames(metadata) <- barcode_ids
        
        # Create Seurat object NOW with filtered data
        incProgress(0.4, detail = "Creating Seurat object with filtered cells...")
        rv$scrna_obj <- CreateSeuratObject(
          counts = count_matrix,
          meta.data = metadata,
          project = "scRNA_analysis"
        )
        
        processing_log <- c(processing_log, paste("Seurat object created with", ncol(rv$scrna_obj), "cells"))
        
        incProgress(0.5, detail = "Filtering genes...")
        seurat_obj <- rv$scrna_obj
        
        # Get initial counts
        n_cells_initial <- ncol(seurat_obj)
        n_genes_initial <- nrow(seurat_obj)
        processing_log <- c(processing_log,
                            "═══ Seurat Object Created ═══",
                            paste("Cells:", n_cells_initial),
                            paste("Genes:", n_genes_initial),
                            "")
        
        incProgress(0.6, detail = "Filtering lowly expressed genes...")
        
        # Step 2: Filter genes
        min_cells <- input$min_cells_per_gene
        
        # Extract counts
        counts <- GetAssayData(seurat_obj, layer = "counts")
        
        # Output a logical vector for every gene on whether there are more than zero counts per cell
        nonzero <- (counts > 0)
        
        # Keep genes expressed in at least min_cells
        keep_genes <- Matrix::rowSums(nonzero) >= min_cells
        n_genes_kept <- sum(keep_genes)
        
        # Get filtered gene names
        filtered_gene_names <- rownames(seurat_obj)[keep_genes]
        filtered_genes_df <- data.frame(
          Index = 1:n_genes_kept,
          Gene_Symbol = filtered_gene_names,
          N_Cells_Expressing = Matrix::rowSums(nonzero)[keep_genes],
          stringsAsFactors = FALSE
        )
        rv$filtered_genes <- filtered_genes_df
        
        processing_log <- c(processing_log,
                            "═══ Gene Filtering ═══",
                            paste("Minimum cells per gene:", min_cells),
                            paste("Genes retained:", n_genes_kept, "/", n_genes_initial,
                                  paste0("(", round(100*n_genes_kept/n_genes_initial, 1), "%)")),
                            paste("Genes removed:", n_genes_initial - n_genes_kept),
                            "")
        
        incProgress(0.6, detail = "Creating filtered Seurat object...")
        
        # Only keep filtered genes
        filtered_counts <- counts[keep_genes, ]
        
        # Create new Seurat object with filtered counts
        seurat_obj_filtered <- CreateSeuratObject(
          counts = filtered_counts,
          meta.data = seurat_obj@meta.data,
          project = "scRNA_filtered"
        )
        
        incProgress(0.7, detail = "Normalizing data...")
        
        # Step 3: Normalize data
        seurat_obj_filtered <- NormalizeData(
          object = seurat_obj_filtered,
          normalization.method = "LogNormalize",
          scale.factor = 10000
        )
        
        processing_log <- c(processing_log,
                            "═══ Normalization ═══",
                            "Method: LogNormalize",
                            "Scale factor: 10000",
                            "✓ Normalization complete",
                            "")
        
        incProgress(0.8, detail = "Scaling data...")
        
        # Scale data
        seurat_obj_filtered <- ScaleData(seurat_obj_filtered, features = rownames(seurat_obj_filtered))
        
        processing_log <- c(processing_log,
                            "✓ Data scaled",
                            "")
        
        incProgress(0.9, detail = "Validating cell type annotations...")
        
        # Step 4: Validate cell type columns
        major_col <- input$celltype_major_column
        minor_col <- input$celltype_minor_column
        
        if (!major_col %in% colnames(seurat_obj_filtered@meta.data)) {
          stop(paste("Major cell type column", major_col, "not found in metadata"))
        }
        
        if (!minor_col %in% colnames(seurat_obj_filtered@meta.data)) {
          stop(paste("Minor cell type column", minor_col, "not found in metadata"))
        }
        
        n_major_types <- length(unique(seurat_obj_filtered@meta.data[[major_col]]))
        n_minor_types <- length(unique(seurat_obj_filtered@meta.data[[minor_col]]))
        
        processing_log <- c(processing_log,
                            "═══ Cell Type Validation ═══",
                            paste("Major cell type column:", major_col),
                            paste("  • Unique types:", n_major_types),
                            paste("Minor cell type column:", minor_col),
                            paste("  • Unique types:", n_minor_types),
                            "")
        
        # Store filtered object
        rv$scrna_obj_filtered <- seurat_obj_filtered
        rv$scrna_processed <- TRUE
        
        incProgress(1, detail = "Complete!")
        
        # Final summary
        processing_log <- c(processing_log,
                            "════════════════════════════════════",
                            "✓✓✓ Preprocessing Complete!",
                            "",
                            "Final Data Summary:",
                            paste("  • Cells:", ncol(seurat_obj_filtered)),
                            paste("  • Genes:", nrow(seurat_obj_filtered)),
                            paste("  • Major cell types:", n_major_types),
                            paste("  • Minor cell types:", n_minor_types),
                            "",
                            "Ready for Step 5: Query Method Analysis",
                            "════════════════════════════════════")
        
        output$scrna_processing_summary <- renderText({
          paste(processing_log, collapse = "\n")
        })
        
      }, error = function(e) {
        output$scrna_processing_summary <- renderText({
          paste("Error during preprocessing:", e$message)
        })
      })
    })
  })
  
  # Run Region-Specific Gene Analysis (Step 2)
  observeEvent(input$run_region_analysis_btn, {
    req(rv$proteomics_validated)
    req(input$group_column)
    
    withProgress(message = 'Analyzing region-specific proteins...', value = 0, {
      
      tryCatch({
        incProgress(0.1, detail = "Preparing data...")
        
        # Extract expression data and metadata
        expr_data <- rv$spatial_expr
        metadata <- rv$metadata
        
        # Ensure first column is gene names
        gene_col <- 1
        gene_names <- expr_data[[gene_col]]
        expr_matrix <- as.matrix(expr_data[, -gene_col])
        rownames(expr_matrix) <- gene_names
        
        # Align samples between expression matrix and metadata
        sample_ids_expr <- colnames(expr_matrix)
        sample_ids_meta <- as.character(metadata[[1]])
        
        # Find common samples and reorder both to match
        common_samples <- intersect(sample_ids_expr, sample_ids_meta)
        
        if (length(common_samples) < 2) {
          stop("Not enough matching samples between expression data and metadata")
        }
        
        # Reorder expression matrix and metadata to match
        expr_matrix <- expr_matrix[, match(common_samples, sample_ids_expr)]
        metadata_ordered <- metadata[match(common_samples, sample_ids_meta), ]
        
        # Get regions/groups
        regions <- unique(metadata_ordered[[input$group_column]])
        regions <- regions[!is.na(regions)]  # Remove NA values
        
        if (length(regions) < 2) {
          stop("Need at least 2 groups for comparison")
        }
        
        incProgress(0.2, detail = "Preparing limma design matrix...")
        
        # Create standard design matrix (one column per group)
        f <- factor(metadata_ordered[[input$group_column]])
        design <- model.matrix(~0 + f)
        colnames(design) <- levels(f)
        
        incProgress(0.3, detail = "Fitting linear model...")
        
        # Fit linear model
        fit <- lmFit(expr_matrix, design)
        
        incProgress(0.4, detail = "Creating contrast matrix...")
        
        # Create contrast matrix for one-vs-rest comparisons
        contrast_list <- list()
        
        for (region in regions) {
          other_regions <- regions[regions != region]
          
          if (length(regions) == 2) {
            # For 2 groups, simple pairwise comparison
            contrast_name <- paste0(region, ".vs.", other_regions[1])
            contrast_formula <- paste0(region, " - ", other_regions[1])
          } else {
            # For >2 groups, compare one vs mean of all others
            contrast_name <- paste0(region, ".vs.Rest")
            # Divide by number of other groups to get proper mean comparison
            n_others <- length(other_regions)
            contrast_formula <- paste0(region, " - (", 
                                       paste(other_regions, collapse = " + "), 
                                       ") / ", n_others)
          }
          
          contrast_list[[contrast_name]] <- contrast_formula
        }
        
        # Build contrast matrix
        contrast.matrix <- makeContrasts(
          contrasts = unlist(contrast_list),
          levels = colnames(design)
        )
        
        incProgress(0.5, detail = "Computing moderated t-statistics...")
        
        # Fit contrasts and compute moderated t-statistics
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit2 <- eBayes(fit2)
        
        incProgress(0.6, detail = "Extracting results...")
        
        # Initialize results list
        region_gene_list <- list()
        all_results <- data.frame()
        
        # Get contrast names (these are the keys from contrast_list)
        contrast_names <- names(contrast_list)
        
        # Extract results for each contrast
        for (i in 1:ncol(coef(fit2))) {
          contrast_name <- contrast_names[i]
          region <- strsplit(contrast_name, "\\.vs\\.")[[1]][1]
          
          # Get top table results
          df <- topTable(fit2, coef = i, adjust.method = "BH", 
                         number = nrow(expr_matrix), sort.by = "P")
          
          # Format results
          region_results <- data.frame(
            Gene = rownames(df),
            Region = region,
            logFC = df$logFC,
            AveExpr = df$AveExpr,
            t = df$t,
            P_value = df$P.Value,
            adj_P_value = df$adj.P.Val,
            Fold_Change = 2^df$logFC,
            stringsAsFactors = FALSE
          )
          
          # Filter significant genes (upregulated in this region)
          # Use adjusted p-value (FDR) instead of raw p-value
          sig_genes <- region_results %>%
            filter(logFC > log2(input$fc_threshold),
                   adj_P_value < input$pval_threshold) %>%
            arrange(adj_P_value)
          
          region_gene_list[[region]] <- sig_genes$Gene
          all_results <- rbind(all_results, sig_genes)
          
          incProgress(0.3 * i / ncol(coef(fit2)), 
                      detail = paste("Processing contrast", i, "of", ncol(coef(fit2))))
        }
        
        incProgress(0.9, detail = "Finalizing results...")
        
        # Store results
        rv$region_specific_genes <- all_results
        rv$region_gene_list <- region_gene_list
        rv$region_analysis_complete <- TRUE
        
        # Generate PCA plot
        incProgress(0.95, detail = "Generating PCA plot...")
        
        # Perform PCA on transposed expression matrix (samples as rows)
        pca_result <- prcomp(t(expr_matrix), scale. = FALSE)
        
        # Create PCA data frame
        pca_data <- data.frame(
          PC1 = pca_result$x[, 1],
          PC2 = pca_result$x[, 2],
          Sample = rownames(pca_result$x),
          Group = metadata_ordered[[input$group_column]]
        )
        
        # Calculate variance explained
        var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
        
        # Store PCA results for download
        rv$pca_result <- pca_result
        rv$pca_data <- pca_data
        rv$var_explained <- var_explained
        
        # Generate PCA plot
        output$pca_plot <- renderPlot({
          # Calculate confidence ellipse for each group
          ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
            geom_point(size = 4, alpha = 0.7) +
            stat_ellipse(level = 0.95, linewidth = 1.2) +
            geom_text_repel(aes(label = Sample), 
                            size = 3.5,
                            box.padding = 0.5,
                            point.padding = 0.3,
                            segment.color = 'grey50',
                            max.overlaps = 20) +
            theme_bw(base_size = 16) +
            labs(
              title = "Principal Component Analysis",
              subtitle = paste0("Samples grouped by ", input$group_column),
              x = paste0("PC1 (", var_explained[1], "% variance)"),
              y = paste0("PC2 (", var_explained[2], "% variance)")
            ) +
            theme(
              legend.position = "right",
              legend.title = element_text(face = "bold"),
              plot.title = element_text(face = "bold", size = 18),
              plot.subtitle = element_text(size = 14)
            ) +
            scale_color_brewer(palette = "Set1")
        })
        
        # Generate Venn diagram for region-specific protein overlaps
        incProgress(0.97, detail = "Generating Venn diagram...")
        
        # Store gene list for Venn diagram
        rv$region_gene_list_for_venn <- region_gene_list
        
        # Calculate overlap statistics
        overlap_stats <- list()
        region_names <- names(region_gene_list)
        
        # Calculate pairwise overlaps
        if (length(region_names) >= 2) {
          overlap_matrix <- matrix(0, nrow = length(region_names), ncol = length(region_names))
          rownames(overlap_matrix) <- region_names
          colnames(overlap_matrix) <- region_names
          
          for (i in 1:length(region_names)) {
            for (j in 1:length(region_names)) {
              if (i == j) {
                overlap_matrix[i, j] <- length(region_gene_list[[i]])
              } else {
                overlap_matrix[i, j] <- length(intersect(region_gene_list[[i]], region_gene_list[[j]]))
              }
            }
          }
          
          rv$overlap_matrix <- overlap_matrix
          
          # Generate overlap statistics text
          stats_text <- paste0("Region-Specific Protein Overlap Statistics\n")
          stats_text <- paste0(stats_text, "==========================================\n\n")
          
          for (region in region_names) {
            stats_text <- paste0(stats_text, region, ": ", length(region_gene_list[[region]]), " proteins\n")
          }
          
          stats_text <- paste0(stats_text, "\nPairwise Overlaps:\n")
          stats_text <- paste0(stats_text, "------------------\n")
          
          for (i in 1:(length(region_names)-1)) {
            for (j in (i+1):length(region_names)) {
              overlap_count <- overlap_matrix[i, j]
              overlap_pct_i <- round(100 * overlap_count / length(region_gene_list[[i]]), 1)
              overlap_pct_j <- round(100 * overlap_count / length(region_gene_list[[j]]), 1)
              stats_text <- paste0(stats_text, 
                                   region_names[i], " ∩ ", region_names[j], ": ",
                                   overlap_count, " proteins ",
                                   "(", overlap_pct_i, "% of ", region_names[i], ", ",
                                   overlap_pct_j, "% of ", region_names[j], ")\n")
            }
          }
          
          # Calculate unique proteins per region
          stats_text <- paste0(stats_text, "\nUnique Proteins (no overlap with other regions):\n")
          stats_text <- paste0(stats_text, "------------------------------------------------\n")
          
          for (region in region_names) {
            all_other_genes <- unique(unlist(region_gene_list[region_names != region]))
            unique_genes <- setdiff(region_gene_list[[region]], all_other_genes)
            unique_pct <- round(100 * length(unique_genes) / length(region_gene_list[[region]]), 1)
            stats_text <- paste0(stats_text, region, ": ", length(unique_genes), 
                                 " proteins (", unique_pct, "%)\n")
          }
          
          stats_text <- paste0(stats_text, "\n💡 Recommendation:\n")
          stats_text <- paste0(stats_text, "------------------\n")
          
          # Calculate average overlap percentage
          total_overlap_pct <- 0
          count_pairs <- 0
          for (i in 1:(length(region_names)-1)) {
            for (j in (i+1):length(region_names)) {
              total_overlap_pct <- total_overlap_pct + 
                (100 * overlap_matrix[i, j] / length(region_gene_list[[i]])) +
                (100 * overlap_matrix[i, j] / length(region_gene_list[[j]]))
              count_pairs <- count_pairs + 2
            }
          }
          avg_overlap_pct <- total_overlap_pct / count_pairs
          
          if (avg_overlap_pct > 30) {
            stats_text <- paste0(stats_text, "⚠️  High overlap detected (average ", round(avg_overlap_pct, 1), "%).\n")
            stats_text <- paste0(stats_text, "   Consider increasing FC threshold or decreasing p-value threshold\n")
            stats_text <- paste0(stats_text, "   to obtain more region-specific markers.\n")
          } else if (avg_overlap_pct > 15) {
            stats_text <- paste0(stats_text, "✓  Moderate overlap (average ", round(avg_overlap_pct, 1), "%).\n")
            stats_text <- paste0(stats_text, "   Current thresholds appear reasonable.\n")
          } else {
            stats_text <- paste0(stats_text, "✓✓ Low overlap (average ", round(avg_overlap_pct, 1), "%).\n")
            stats_text <- paste0(stats_text, "   Excellent region specificity!\n")
          }
          
          output$overlap_stats <- renderText({ stats_text })
        }
        
        # Generate Venn diagram plot
        output$venn_plot <- renderPlot({
          # Limit to first 5 regions if more (VennDiagram supports up to 5)
          if (length(region_gene_list) > 5) {
            showNotification("Note: Venn diagram limited to first 5 regions. All regions included in statistics.", 
                             type = "warning", duration = 5)
            plot_list <- region_gene_list[1:5]
          } else {
            plot_list <- region_gene_list
          }
          
          if (length(plot_list) < 2) {
            plot.new()
            text(0.5, 0.5, "Need at least 2 regions for Venn diagram", cex = 1.5)
          } else {
            # Try to create Venn diagram based on number of sets
            n_sets <- length(plot_list)
            
            tryCatch({
              # Use VennDiagram package
              colors <- c("#FF6B6B", "#4ECDC4", "#45B7D1", "#FFA07A", "#98D8C8")[1:n_sets]
              
              if (n_sets == 2) {
                venn_plot <- venn.diagram(
                  x = plot_list,
                  category.names = names(plot_list),
                  filename = NULL,
                  output = TRUE,
                  imagetype = "png",
                  height = 2000,
                  width = 2000,
                  resolution = 300,
                  compression = "lzw",
                  lwd = 2,
                  col = "black",
                  fill = colors,
                  cex = 1.5,
                  fontface = "bold",
                  fontfamily = "sans",
                  cat.cex = 1.5,
                  cat.fontface = "bold",
                  cat.default.pos = "outer",
                  cat.fontfamily = "sans",
                  main = "Overlap of Region-Specific Proteins",
                  main.cex = 1.8,
                  sub = paste0("FC > ", input$fc_threshold, ", Adj.P < ", input$pval_threshold),
                  sub.cex = 1.2
                )
              } else if (n_sets == 3) {
                venn_plot <- venn.diagram(
                  x = plot_list,
                  category.names = names(plot_list),
                  filename = NULL,
                  output = TRUE,
                  imagetype = "png",
                  height = 2000,
                  width = 2000,
                  resolution = 300,
                  compression = "lzw",
                  lwd = 2,
                  col = "black",
                  fill = colors,
                  cex = 1.5,
                  fontface = "bold",
                  fontfamily = "sans",
                  cat.cex = 1.5,
                  cat.fontface = "bold",
                  cat.fontfamily = "sans",
                  main = "Overlap of Region-Specific Proteins",
                  main.cex = 1.8,
                  sub = paste0("FC > ", input$fc_threshold, ", Adj.P < ", input$pval_threshold),
                  sub.cex = 1.2
                )
              } else if (n_sets == 4) {
                venn_plot <- venn.diagram(
                  x = plot_list,
                  category.names = names(plot_list),
                  filename = NULL,
                  output = TRUE,
                  imagetype = "png",
                  height = 2000,
                  width = 2000,
                  resolution = 300,
                  compression = "lzw",
                  lwd = 2,
                  col = "black",
                  fill = colors,
                  cex = 1.2,
                  fontface = "bold",
                  fontfamily = "sans",
                  cat.cex = 1.2,
                  cat.fontface = "bold",
                  cat.fontfamily = "sans",
                  main = "Overlap of Region-Specific Proteins",
                  main.cex = 1.8,
                  sub = paste0("FC > ", input$fc_threshold, ", Adj.P < ", input$pval_threshold),
                  sub.cex = 1.2
                )
              } else if (n_sets == 5) {
                venn_plot <- venn.diagram(
                  x = plot_list,
                  category.names = names(plot_list),
                  filename = NULL,
                  output = TRUE,
                  imagetype = "png",
                  height = 2000,
                  width = 2000,
                  resolution = 300,
                  compression = "lzw",
                  lwd = 2,
                  col = "black",
                  fill = colors,
                  cex = 1,
                  fontface = "bold",
                  fontfamily = "sans",
                  cat.cex = 1,
                  cat.fontface = "bold",
                  cat.fontfamily = "sans",
                  main = "Overlap of Region-Specific Proteins",
                  main.cex = 1.8,
                  sub = paste0("FC > ", input$fc_threshold, ", Adj.P < ", input$pval_threshold),
                  sub.cex = 1.2
                )
              }
              
              grid.newpage()
              grid.draw(venn_plot)
              
            }, error = function(e) {
              # Fallback: Create a simple ggplot2 visualization showing overlaps
              overlap_data <- data.frame(
                Region = character(),
                Unique = numeric(),
                Shared = numeric(),
                stringsAsFactors = FALSE
              )
              
              for (region in names(plot_list)) {
                all_other_genes <- unique(unlist(plot_list[names(plot_list) != region]))
                unique_genes <- length(setdiff(plot_list[[region]], all_other_genes))
                shared_genes <- length(intersect(plot_list[[region]], all_other_genes))
                
                overlap_data <- rbind(overlap_data, data.frame(
                  Region = region,
                  Unique = unique_genes,
                  Shared = shared_genes,
                  stringsAsFactors = FALSE
                ))
              }
              
              # Reshape for stacked bar plot
              overlap_long <- reshape2::melt(overlap_data, id.vars = "Region", 
                                             variable.name = "Type", value.name = "Count")
              
              ggplot(overlap_long, aes(x = Region, y = Count, fill = Type)) +
                geom_bar(stat = "identity") +
                scale_fill_manual(values = c("Unique" = "#4ECDC4", "Shared" = "#FF6B6B"),
                                  labels = c("Unique to Region", "Shared with Others")) +
                theme_minimal(base_size = 14) +
                labs(title = "Overlap of Region-Specific Proteins",
                     subtitle = paste0("FC > ", input$fc_threshold, ", Adj.P < ", input$pval_threshold, 
                                       "\n(Venn diagram unavailable - showing bar plot)"),
                     x = "Region",
                     y = "Number of Proteins",
                     fill = "") +
                theme(axis.text.x = element_text(angle = 45, hjust = 1),
                      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
                      plot.subtitle = element_text(size = 12, hjust = 0.5),
                      legend.position = "top")
            })
          }
        })
        
        incProgress(1, detail = "Complete!")
        
        output$region_analysis_status <- renderUI({
          div(
            style = "background-color: #dff0d8; padding: 15px; border-radius: 5px; margin: 10px 0;",
            h4(icon("check-circle"), "Analysis Complete!", style = "color: green;"),
            p(paste("Identified", nrow(rv$region_specific_genes), "region-specific proteins across", 
                    length(region_gene_list), "regions")),
            p("Proteins per region:"),
            tags$ul(
              lapply(names(region_gene_list), function(r) {
                tags$li(paste(r, "> all others:", length(region_gene_list[[r]]), "proteins"))
              })
            )
          )
        })
        
        output$region_genes_table <- renderDT({
          datatable(rv$region_specific_genes %>% 
                      mutate(logFC = round(logFC, 3),
                             AveExpr = round(AveExpr, 2),
                             t = round(t, 3),
                             P_value = formatC(P_value, format = "e", digits = 2),
                             adj_P_value = formatC(adj_P_value, format = "e", digits = 2),
                             Fold_Change = round(Fold_Change, 2)), 
                    options = list(pageLength = 10, scrollX = TRUE),
                    filter = 'top')
        })
        
        showNotification("Region-specific protein analysis complete (limma moderated t-test)!", type = "message", duration = 5)
        
      }, error = function(e) {
        output$region_analysis_status <- renderUI({
          div(
            style = "background-color: #f2dede; padding: 15px; border-radius: 5px;",
            h4(icon("exclamation-circle"), "Error:", style = "color: red;"),
            p(e$message)
          )
        })
        showNotification(paste("Error:", e$message), type = "error", duration = 10)
      })
    })
  })
  
  # Run Query Method Analysis (Step 5)
  observeEvent(input$run_query_method_btn, {
    req(rv$region_analysis_complete)
    req(rv$scrna_processed)
    req(rv$scrna_obj_filtered)
    
    withProgress(message = 'Running Query Method Analysis...', value = 0, {
      
      tryCatch({
        incProgress(0.1, detail = "Preparing scRNA-seq data...")
        
        # Get scRNA-seq data (use filtered/normalized object)
        scrna_obj <- rv$scrna_obj_filtered
        
        # Ensure cell type annotation exists - use the major/minor columns from Step 4
        minor_col <- input$celltype_minor_column
        major_col <- input$celltype_major_column
        if (is.null(minor_col) || !minor_col %in% colnames(scrna_obj@meta.data)) {
          stop("Minor cell type column not found. Please complete Step 4 preprocessing first.")
        }
        if (is.null(major_col) || !major_col %in% colnames(scrna_obj@meta.data)) {
          stop("Major cell type column not found. Please complete Step 4 preprocessing first.")
        }
        
        scrna_obj@meta.data$celltype <- scrna_obj@meta.data[[minor_col]]
        
        # Get cell types
        cell_types <- unique(scrna_obj@meta.data$celltype)
        
        # Get gene expression matrix (normalized data)
        expr_matrix <- GetAssayData(scrna_obj, layer = "data")
        
        incProgress(0.2, detail = "Matching genes...")
        
        # Build overlap between proteomics DEPs and scRNA-seq genes
        scrna_genes <- if (!is.null(rv$filtered_genes) && "Gene_Symbol" %in% colnames(rv$filtered_genes)) {
          rv$filtered_genes$Gene_Symbol
        } else {
          rownames(expr_matrix)
        }
        
        rv$region_gene_list_overlap <- lapply(rv$region_gene_list, function(genes) {
          intersect(genes, scrna_genes)
        })
        
        if (length(rv$region_gene_list) == 0) {
          overlap_summary <- data.frame(
            Region = character(),
            Proteomics_DEPs = integer(),
            Overlap_scRNAseq = integer(),
            stringsAsFactors = FALSE
          )
        } else {
          region_names <- names(rv$region_gene_list)
          if (is.null(region_names) || any(region_names == "")) {
            region_names <- paste0("Region_", seq_along(rv$region_gene_list))
          }
          names(rv$region_gene_list_overlap) <- region_names
          overlap_summary <- data.frame(
            Region = region_names,
            Proteomics_DEPs = vapply(rv$region_gene_list, length, integer(1)),
            Overlap_scRNAseq = vapply(rv$region_gene_list_overlap, length, integer(1)),
            stringsAsFactors = FALSE
          )
        }
        
        overlap_genes_list <- lapply(names(rv$region_gene_list_overlap), function(region) {
          genes <- rv$region_gene_list_overlap[[region]]
          if (length(genes) == 0) return(NULL)
          data.frame(Region = region, Gene = genes, stringsAsFactors = FALSE)
        })
        
        overlap_genes_list <- overlap_genes_list[!sapply(overlap_genes_list, is.null)]
        
        if (length(overlap_genes_list) == 0) {
          overlap_genes_df <- data.frame(Region = character(), Gene = character(), stringsAsFactors = FALSE)
        } else {
          overlap_genes_df <- do.call(rbind, overlap_genes_list)
        }
        
        output$overlap_summary_table <- renderDT({
          datatable(overlap_summary, options = list(pageLength = 10, scrollX = TRUE))
        })
        
        output$overlap_genes_table <- renderDT({
          datatable(overlap_genes_df, options = list(pageLength = 10, scrollX = TRUE), filter = 'top')
        })
        
        # Initialize results
        query_results_list <- list()
        
        # For each region
        regions <- names(rv$region_gene_list)
        
        for (region_idx in seq_along(regions)) {
          region <- regions[region_idx]
          region_genes <- rv$region_gene_list_overlap[[region]]
          
          # Match genes with scRNA-seq data
          matched_genes <- intersect(region_genes, rownames(expr_matrix))
          
          if (length(matched_genes) == 0) {
            warning(paste("No matching genes for region", region))
            next
          }
          
          incProgress(0.3 * region_idx / length(regions), 
                      detail = paste("Analyzing region", region))
          
          # Extract expression for matched genes
          region_expr <- expr_matrix[matched_genes, , drop = FALSE]
          
          # Normalize expression per gene (z-score)
          region_expr_norm <- t(scale(t(as.matrix(region_expr))))
          
          # Calculate enrichment score for each cell type
          enrichment_scores <- data.frame(
            Region = character(),
            Major_Cell_Type = character(),
            Minor_Cell_Type = character(),
            Signature_Score = numeric(),
            Score_Per_Gene = numeric(),
            N_Cells = integer(),
            N_Genes = integer(),
            stringsAsFactors = FALSE
          )
          
          for (ct in cell_types) {
            # Get cells of this type
            cells_ct <- colnames(scrna_obj)[scrna_obj@meta.data$celltype == ct]
            
            if (length(cells_ct) == 0) next
            
            # Extract z-normalized expression for this cell type
            cell_type_data <- region_expr_norm[, cells_ct, drop = FALSE]
            
            # Signature score = sum(z_scores) / n_cells (matches user's method)
            signature_score <- sum(cell_type_data, na.rm = TRUE) / length(cells_ct)
            
            # Score per gene
            score_per_gene <- signature_score / length(matched_genes)
            
            major_ct <- scrna_obj@meta.data[cells_ct, major_col]
            major_ct <- major_ct[!is.na(major_ct)]
            major_ct <- if (length(major_ct) == 0) NA else major_ct[1]
            
            enrichment_scores <- rbind(enrichment_scores, data.frame(
              Region = region,
              Major_Cell_Type = as.character(major_ct),
              Minor_Cell_Type = ct,
              Signature_Score = signature_score,
              Score_Per_Gene = score_per_gene,
              N_Cells = length(cells_ct),
              N_Genes = length(matched_genes),
              stringsAsFactors = FALSE
            ))
          }
          
          query_results_list[[region]] <- enrichment_scores
        }
        
        incProgress(0.8, detail = "Combining results...")
        
        # Combine all results
        if (length(query_results_list) == 0) {
          rv$query_results <- data.frame(
            Region = character(),
            Major_Cell_Type = character(),
            Minor_Cell_Type = character(),
            Signature_Score = numeric(),
            Score_Per_Gene = numeric(),
            N_Cells = integer(),
            N_Genes = integer(),
            stringsAsFactors = FALSE
          )
        } else {
          rv$query_results <- do.call(rbind, query_results_list)
        }
        
        # Rank by signature score for easier interpretation
        if (nrow(rv$query_results) > 0) {
          rv$query_results <- rv$query_results[order(-rv$query_results$Signature_Score), ]
        }
        
        incProgress(0.9, detail = "Creating visualizations...")
        
        rv$query_method_complete <- TRUE
        
        # Build cell type ordering (major -> minor, alphabetical)
        minor_col <- input$celltype_minor_column
        major_col <- input$celltype_major_column
        celltype_map <- unique(scrna_obj@meta.data[, c(major_col, minor_col), drop = FALSE])
        celltype_map <- celltype_map[order(celltype_map[[major_col]], celltype_map[[minor_col]]), , drop = FALSE]
        celltype_levels <- celltype_map[[minor_col]]
        
        # Create positive enrichment plot
        output$query_bubble_plot_positive <- renderPlot({
          plot_data <- rv$query_results[rv$query_results$Signature_Score > 0, , drop = FALSE]
          plot_data$Minor_Cell_Type <- factor(plot_data$Minor_Cell_Type, levels = celltype_levels)
          
          if (nrow(plot_data) == 0) {
            plot.new()
            text(0.5, 0.5, "No positive enrichment", cex = 1.5)
            return()
          }
          
          ggplot(plot_data, aes(x = Region, y = Minor_Cell_Type)) +
            geom_point(aes(size = abs(Signature_Score), color = Signature_Score)) +
            scale_color_gradient2(low = "white", mid = "white", high = "red",
                                  midpoint = 0, name = "Signature Score") +
            scale_size_continuous(range = c(2, 15), name = "|Score|") +
            theme_minimal(base_size = 14) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                  axis.text.y = element_text(size = 12),
                  panel.grid.major = element_line(color = "gray90"),
                  panel.border = element_rect(fill = NA, color = "black")) +
            labs(title = "Positive Signature Scores",
                 x = "Histopathological Region",
                 y = "Cell Type")
        })
        
        # Create negative enrichment plot
        output$query_bubble_plot_negative <- renderPlot({
          plot_data <- rv$query_results[rv$query_results$Signature_Score < 0, , drop = FALSE]
          plot_data$Minor_Cell_Type <- factor(plot_data$Minor_Cell_Type, levels = celltype_levels)
          
          if (nrow(plot_data) == 0) {
            plot.new()
            text(0.5, 0.5, "No negative enrichment", cex = 1.5)
            return()
          }
          
          ggplot(plot_data, aes(x = Region, y = Minor_Cell_Type)) +
            geom_point(aes(size = abs(Signature_Score), color = Signature_Score)) +
            scale_color_gradient2(low = "blue", mid = "white", high = "white",
                                  midpoint = 0, name = "Signature Score") +
            scale_size_continuous(range = c(2, 15), name = "|Score|") +
            theme_minimal(base_size = 14) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                  axis.text.y = element_text(size = 12),
                  panel.grid.major = element_line(color = "gray90"),
                  panel.border = element_rect(fill = NA, color = "black")) +
            labs(title = "Negative Signature Scores",
                 x = "Histopathological Region",
                 y = "Cell Type")
        })
        
        output$query_results_table <- renderDT({
          datatable(rv$query_results %>%
                      mutate(Signature_Score = round(Signature_Score, 4),
                             Score_Per_Gene = round(Score_Per_Gene, 6)),
                    options = list(pageLength = 20, scrollX = TRUE),
                    filter = 'top',
                    rownames = FALSE)
        })
        
        output$query_method_status <- renderUI({
          top_positive <- rv$query_results[rv$query_results$Signature_Score > 0, ]
          top_negative <- rv$query_results[rv$query_results$Signature_Score < 0, ]
          
          div(
            style = "background-color: #dff0d8; padding: 15px; border-radius: 5px; margin: 10px 0;",
            h4(icon("check-circle"), "Signature Score Analysis Complete!", style = "color: green;"),
            p(paste("Total region-cell type combinations:", nrow(rv$query_results))),
            p(paste("Positive scores (enriched):", nrow(top_positive))),
            p(paste("Negative scores (depleted):", nrow(top_negative))),
            p("Scroll down to view bubble plot and detailed results.")
          )
        })
        
        incProgress(1, detail = "Complete!")
        
        showNotification("Query Method analysis complete!", type = "message", duration = 5)
        
      }, error = function(e) {
        output$query_method_status <- renderUI({
          div(
            style = "background-color: #f2dede; padding: 15px; border-radius: 5px;",
            h4(icon("exclamation-circle"), "Error:", style = "color: red;"),
            p(e$message)
          )
        })
        showNotification(paste("Error:", e$message), type = "error", duration = 10)
      })
    })
  })
  
  # Download handlers
  output$download_region_genes <- downloadHandler(
    filename = function() {
      paste0("region_specific_genes_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(rv$region_specific_genes, file, row.names = FALSE)
    }
  )
  
  output$download_overlap_genes <- downloadHandler(
    filename = function() {
      paste0("overlap_genes_", Sys.Date(), ".csv")
    },
    content = function(file) {
      if (is.null(rv$region_gene_list_overlap)) {
        write.csv(data.frame(), file, row.names = FALSE)
      } else {
        overlap_genes_df <- do.call(rbind, lapply(names(rv$region_gene_list_overlap), function(region) {
          genes <- rv$region_gene_list_overlap[[region]]
          if (length(genes) == 0) return(NULL)
          data.frame(Region = region, Gene = genes, stringsAsFactors = FALSE)
        }))
        if (is.null(overlap_genes_df)) overlap_genes_df <- data.frame(Region = character(), Gene = character())
        write.csv(overlap_genes_df, file, row.names = FALSE)
      }
    }
  )
  
  output$download_pca_plot <- downloadHandler(
    filename = function() {
      paste0("pca_plot_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      pdf(file, width = 10, height = 8)
      
      # Recreate the PCA plot for download
      p <- ggplot(rv$pca_data, aes(x = PC1, y = PC2, color = Group)) +
        geom_point(size = 4, alpha = 0.7) +
        stat_ellipse(level = 0.95, linewidth = 1.2) +
        geom_text_repel(aes(label = Sample), 
                        size = 3.5,
                        box.padding = 0.5,
                        point.padding = 0.3,
                        segment.color = 'grey50',
                        max.overlaps = 20) +
        theme_bw(base_size = 16) +
        labs(
          title = "Principal Component Analysis",
          subtitle = paste0("Samples grouped by ", input$group_column),
          x = paste0("PC1 (", rv$var_explained[1], "% variance)"),
          y = paste0("PC2 (", rv$var_explained[2], "% variance)")
        ) +
        theme(
          legend.position = "right",
          legend.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = 18),
          plot.subtitle = element_text(size = 14)
        ) +
        scale_color_brewer(palette = "Set1")
      
      print(p)
      dev.off()
    }
  )
  
  output$download_venn_plot <- downloadHandler(
    filename = function() {
      paste0("venn_diagram_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      # Recreate Venn diagram for download
      plot_list <- rv$region_gene_list_for_venn
      
      if (length(plot_list) > 5) {
        plot_list <- plot_list[1:5]
      }
      
      if (length(plot_list) >= 2) {
        n_sets <- length(plot_list)
        colors <- c("#FF6B6B", "#4ECDC4", "#45B7D1", "#FFA07A", "#98D8C8")[1:n_sets]
        
        tryCatch({
          # Use VennDiagram package
          if (n_sets == 2) {
            venn_plot <- venn.diagram(
              x = plot_list,
              category.names = names(plot_list),
              filename = file,
              imagetype = "pdf",
              height = 8,
              width = 10,
              lwd = 2,
              col = "black",
              fill = colors,
              cex = 1.5,
              fontface = "bold",
              fontfamily = "sans",
              cat.cex = 1.5,
              cat.fontface = "bold",
              cat.default.pos = "outer",
              cat.fontfamily = "sans",
              main = "Overlap of Region-Specific Proteins",
              main.cex = 1.8,
              sub = paste0("FC > ", input$fc_threshold, ", Adj.P < ", input$pval_threshold),
              sub.cex = 1.2
            )
          } else if (n_sets == 3) {
            venn_plot <- venn.diagram(
              x = plot_list,
              category.names = names(plot_list),
              filename = file,
              imagetype = "pdf",
              height = 8,
              width = 10,
              lwd = 2,
              col = "black",
              fill = colors,
              cex = 1.5,
              fontface = "bold",
              fontfamily = "sans",
              cat.cex = 1.5,
              cat.fontface = "bold",
              cat.fontfamily = "sans",
              main = "Overlap of Region-Specific Proteins",
              main.cex = 1.8,
              sub = paste0("FC > ", input$fc_threshold, ", Adj.P < ", input$pval_threshold),
              sub.cex = 1.2
            )
          } else if (n_sets == 4) {
            venn_plot <- venn.diagram(
              x = plot_list,
              category.names = names(plot_list),
              filename = file,
              imagetype = "pdf",
              height = 8,
              width = 10,
              lwd = 2,
              col = "black",
              fill = colors,
              cex = 1.2,
              fontface = "bold",
              fontfamily = "sans",
              cat.cex = 1.2,
              cat.fontface = "bold",
              cat.fontfamily = "sans",
              main = "Overlap of Region-Specific Proteins",
              main.cex = 1.8,
              sub = paste0("FC > ", input$fc_threshold, ", Adj.P < ", input$pval_threshold),
              sub.cex = 1.2
            )
          } else if (n_sets == 5) {
            venn_plot <- venn.diagram(
              x = plot_list,
              category.names = names(plot_list),
              filename = file,
              imagetype = "pdf",
              height = 8,
              width = 10,
              lwd = 2,
              col = "black",
              fill = colors,
              cex = 1,
              fontface = "bold",
              fontfamily = "sans",
              cat.cex = 1,
              cat.fontface = "bold",
              cat.fontfamily = "sans",
              main = "Overlap of Region-Specific Proteins",
              main.cex = 1.8,
              sub = paste0("FC > ", input$fc_threshold, ", Adj.P < ", input$pval_threshold),
              sub.cex = 1.2
            )
          }
        }, error = function(e) {
          # Fallback: Create a simple bar plot
          pdf(file, width = 10, height = 8)
          
          overlap_data <- data.frame(
            Region = character(),
            Unique = numeric(),
            Shared = numeric(),
            stringsAsFactors = FALSE
          )
          
          for (region in names(plot_list)) {
            all_other_genes <- unique(unlist(plot_list[names(plot_list) != region]))
            unique_genes <- length(setdiff(plot_list[[region]], all_other_genes))
            shared_genes <- length(intersect(plot_list[[region]], all_other_genes))
            
            overlap_data <- rbind(overlap_data, data.frame(
              Region = region,
              Unique = unique_genes,
              Shared = shared_genes,
              stringsAsFactors = FALSE
            ))
          }
          
          overlap_long <- reshape2::melt(overlap_data, id.vars = "Region", 
                                         variable.name = "Type", value.name = "Count")
          
          p <- ggplot(overlap_long, aes(x = Region, y = Count, fill = Type)) +
            geom_bar(stat = "identity") +
            scale_fill_manual(values = c("Unique" = "#4ECDC4", "Shared" = "#FF6B6B"),
                              labels = c("Unique to Region", "Shared with Others")) +
            theme_minimal(base_size = 14) +
            labs(title = "Overlap of Region-Specific Proteins",
                 subtitle = paste0("FC > ", input$fc_threshold, ", Adj.P < ", input$pval_threshold),
                 x = "Region",
                 y = "Number of Proteins",
                 fill = "") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
                  plot.subtitle = element_text(size = 12, hjust = 0.5),
                  legend.position = "top")
          
          print(p)
          dev.off()
        })
      }
    }
  )
  
  output$download_query_results <- downloadHandler(
    filename = function() {
      paste0("query_method_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(rv$query_results, file, row.names = FALSE)
    }
  )
  
  output$download_query_plot <- downloadHandler(
    filename = function() {
      paste0("query_bubble_plot_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      pdf(file, width = 12, height = 10)
      plot_data_pos <- rv$query_results[rv$query_results$Signature_Score > 0, , drop = FALSE]
      plot_data_neg <- rv$query_results[rv$query_results$Signature_Score < 0, , drop = FALSE]
      
      if (nrow(plot_data_pos) > 0) {
        p_pos <- ggplot(plot_data_pos, aes(x = Region, y = Minor_Cell_Type)) +
          geom_point(aes(size = abs(Signature_Score), color = Signature_Score)) +
          scale_color_gradient2(low = "white", mid = "white", high = "red",
                                midpoint = 0, name = "Signature Score") +
          scale_size_continuous(range = c(2, 15), name = "|Score|") +
          theme_minimal(base_size = 14) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                axis.text.y = element_text(size = 12),
                panel.grid.major = element_line(color = "gray90"),
                panel.border = element_rect(fill = NA, color = "black")) +
          labs(title = "Positive Signature Scores",
               x = "Histopathological Region",
               y = "Cell Type")
        print(p_pos)
      } else {
        plot.new()
        text(0.5, 0.5, "No positive enrichment", cex = 1.5)
      }
      
      if (nrow(plot_data_neg) > 0) {
        p_neg <- ggplot(plot_data_neg, aes(x = Region, y = Minor_Cell_Type)) +
          geom_point(aes(size = abs(Signature_Score), color = Signature_Score)) +
          scale_color_gradient2(low = "blue", mid = "white", high = "white",
                                midpoint = 0, name = "Signature Score") +
          scale_size_continuous(range = c(2, 15), name = "|Score|") +
          theme_minimal(base_size = 14) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                axis.text.y = element_text(size = 12),
                panel.grid.major = element_line(color = "gray90"),
                panel.border = element_rect(fill = NA, color = "black")) +
          labs(title = "Negative Signature Scores",
               x = "Histopathological Region",
               y = "Cell Type")
        print(p_neg)
      } else {
        plot.new()
        text(0.5, 0.5, "No negative enrichment", cex = 1.5)
      }
      
      dev.off()
    }
  )

  output$download_all_results <- downloadHandler(
    filename = function() {
      paste0("all_results_", Sys.Date(), ".zip")
    },
    content = function(file) {
      tmp_dir <- tempdir()
      files_to_zip <- c()
      
      # Save region-specific genes table
      if (!is.null(rv$region_specific_genes)) {
        f <- file.path(tmp_dir, "region_specific_genes.csv")
        write.csv(rv$region_specific_genes, f, row.names = FALSE)
        files_to_zip <- c(files_to_zip, "region_specific_genes.csv")
      }
      
      # Save overlap summary table
      if (!is.null(rv$region_gene_list) && !is.null(rv$region_gene_list_overlap)) {
        region_names <- names(rv$region_gene_list)
        if (is.null(region_names) || any(region_names == "")) {
          region_names <- paste0("Region_", seq_along(rv$region_gene_list))
        }
        overlap_summary <- data.frame(
          Region = region_names,
          Proteomics_DEPs = vapply(rv$region_gene_list, length, integer(1)),
          Overlap_scRNAseq = vapply(rv$region_gene_list_overlap, length, integer(1)),
          stringsAsFactors = FALSE
        )
        f <- file.path(tmp_dir, "overlap_summary.csv")
        write.csv(overlap_summary, f, row.names = FALSE)
        files_to_zip <- c(files_to_zip, "overlap_summary.csv")
      }
      
      # Save overlap genes table
      if (!is.null(rv$region_gene_list_overlap)) {
        overlap_genes_df <- do.call(rbind, lapply(names(rv$region_gene_list_overlap), function(region) {
          genes <- rv$region_gene_list_overlap[[region]]
          if (length(genes) == 0) return(NULL)
          data.frame(Region = region, Gene = genes, stringsAsFactors = FALSE)
        }))
        if (!is.null(overlap_genes_df)) {
          f <- file.path(tmp_dir, "overlap_genes.csv")
          write.csv(overlap_genes_df, f, row.names = FALSE)
          files_to_zip <- c(files_to_zip, "overlap_genes.csv")
        }
      }
      
      # Save query method results table
      if (!is.null(rv$query_results)) {
        f <- file.path(tmp_dir, "query_method_results.csv")
        write.csv(rv$query_results, f, row.names = FALSE)
        files_to_zip <- c(files_to_zip, "query_method_results.csv")
      }
      
      if (length(files_to_zip) == 0) {
        f <- file.path(tmp_dir, "empty_results.txt")
        writeLines("No results available yet.", f)
        files_to_zip <- c(files_to_zip, "empty_results.txt")
      }
      
      # Change to temp directory and zip files without full paths
      old_wd <- getwd()
      setwd(tmp_dir)
      utils::zip(zipfile = file, files = files_to_zip)
      setwd(old_wd)
    }
  )
  
  output$download_report <- downloadHandler(
    filename = function() {
      paste0("summary_report_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      pdf(file, width = 11, height = 8.5)
      
      # Page 1: Summary information
      plot.new()
      text(0.5, 0.95, "Multi-Omics Integration Analysis Report", adj = 0.5, cex = 1.8, font = 2)
      text(0.5, 0.90, paste("Date:", Sys.Date()), adj = 0.5, cex = 1.2)
      text(0.1, 0.80, "Analysis Summary:", adj = 0, cex = 1.4, font = 2)
      y_pos <- 0.75
      if (rv$proteomics_validated) {
        text(0.1, y_pos, "[OK] Step 1: Proteomics Data Validated", adj = 0, cex = 1.1)
        y_pos <- y_pos - 0.05
      }
      if (!is.null(rv$region_specific_genes)) {
        text(0.1, y_pos, paste0("[OK] Step 2: ", nrow(rv$region_specific_genes), " Region-Specific Proteins Identified"), adj = 0, cex = 1.1)
        y_pos <- y_pos - 0.05
      }
      if (rv$scrna_loaded) {
        text(0.1, y_pos, "[OK] Step 3: scRNA-seq Data Loaded", adj = 0, cex = 1.1)
        y_pos <- y_pos - 0.05
      }
      if (rv$scrna_processed) {
        text(0.1, y_pos, "[OK] Step 4: scRNA-seq Data Processed", adj = 0, cex = 1.1)
        y_pos <- y_pos - 0.05
      }
      if (!is.null(rv$query_results)) {
        text(0.1, y_pos, paste0("[OK] Step 5: ", nrow(rv$query_results), " Query Method Results Generated"), adj = 0, cex = 1.1)
        y_pos <- y_pos - 0.05
      }
      
      # Page 2: PCA Plot
      if (!is.null(rv$pca_data) && !is.null(rv$var_explained)) {
        p <- ggplot(rv$pca_data, aes(x = PC1, y = PC2, color = Group)) +
          geom_point(size = 4, alpha = 0.7) +
          stat_ellipse(level = 0.95, linewidth = 1.2) +
          geom_text_repel(aes(label = Sample), 
                          size = 3.5,
                          box.padding = 0.5,
                          point.padding = 0.3,
                          segment.color = 'grey50',
                          max.overlaps = 20) +
          theme_bw(base_size = 16) +
          labs(
            title = "Principal Component Analysis",
            x = paste0("PC1 (", rv$var_explained[1], "% variance)"),
            y = paste0("PC2 (", rv$var_explained[2], "% variance)")
          ) +
          theme(
            legend.position = "right",
            legend.title = element_text(face = "bold"),
            plot.title = element_text(face = "bold", size = 18, hjust = 0.5)
          ) +
          scale_color_brewer(palette = "Set1")
        print(p)
      }
      
      # Page 3: Venn Diagram
      if (!is.null(rv$region_gene_list_for_venn)) {
        plot_list <- rv$region_gene_list_for_venn
        if (length(plot_list) > 5) plot_list <- plot_list[1:5]
        
        if (length(plot_list) >= 2) {
          n_sets <- length(plot_list)
          colors <- c("#FF6B6B", "#4ECDC4", "#45B7D1", "#FFA07A", "#98D8C8")[1:n_sets]
          
          tryCatch({
            grid.newpage()
            if (n_sets == 2) {
              venn_plot <- venn.diagram(
                x = plot_list,
                category.names = names(plot_list),
                filename = NULL,
                lwd = 2,
                col = "black",
                fill = colors,
                cex = 1.5,
                fontface = "bold",
                cat.cex = 1.5,
                cat.fontface = "bold",
                cat.default.pos = "outer",
                main = "Overlap of Region-Specific Proteins",
                main.cex = 1.5
              )
            } else if (n_sets == 3) {
              venn_plot <- venn.diagram(
                x = plot_list,
                category.names = names(plot_list),
                filename = NULL,
                lwd = 2,
                col = "black",
                fill = colors,
                cex = 1.5,
                fontface = "bold",
                cat.cex = 1.5,
                cat.fontface = "bold",
                main = "Overlap of Region-Specific Proteins",
                main.cex = 1.5
              )
            } else if (n_sets >= 4) {
              venn_plot <- venn.diagram(
                x = plot_list,
                category.names = names(plot_list),
                filename = NULL,
                lwd = 2,
                col = "black",
                fill = colors,
                cex = 1,
                fontface = "bold",
                cat.cex = 1,
                cat.fontface = "bold",
                main = "Overlap of Region-Specific Proteins",
                main.cex = 1.5
              )
            }
            grid.draw(venn_plot)
          }, error = function(e) {
            # Fallback bar plot
            overlap_data <- data.frame(
              Region = character(),
              Unique = numeric(),
              Shared = numeric(),
              stringsAsFactors = FALSE
            )
            
            for (region in names(plot_list)) {
              all_other_genes <- unique(unlist(plot_list[names(plot_list) != region]))
              unique_genes <- length(setdiff(plot_list[[region]], all_other_genes))
              shared_genes <- length(intersect(plot_list[[region]], all_other_genes))
              
              overlap_data <- rbind(overlap_data, data.frame(
                Region = region,
                Unique = unique_genes,
                Shared = shared_genes,
                stringsAsFactors = FALSE
              ))
            }
            
            overlap_long <- reshape2::melt(overlap_data, id.vars = "Region", 
                                           variable.name = "Type", value.name = "Count")
            
            p <- ggplot(overlap_long, aes(x = Region, y = Count, fill = Type)) +
              geom_bar(stat = "identity") +
              scale_fill_manual(values = c("Unique" = "#4ECDC4", "Shared" = "#FF6B6B"),
                                labels = c("Unique to Region", "Shared with Others")) +
              theme_minimal(base_size = 14) +
              labs(title = "Overlap of Region-Specific Proteins",
                   x = "Region",
                   y = "Number of Proteins",
                   fill = "") +
              theme(axis.text.x = element_text(angle = 45, hjust = 1),
                    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
                    legend.position = "top")
            print(p)
          })
        }
      }
      
      # Page 4: Positive Enrichment Bubble Plot
      if (!is.null(rv$query_results)) {
        plot_data_pos <- rv$query_results[rv$query_results$Signature_Score > 0, , drop = FALSE]
        
        if (nrow(plot_data_pos) > 0) {
          p_pos <- ggplot(plot_data_pos, aes(x = Region, y = Minor_Cell_Type)) +
            geom_point(aes(size = abs(Signature_Score), color = Signature_Score)) +
            scale_color_gradient2(low = "white", mid = "white", high = "red",
                                  midpoint = 0, name = "Signature Score") +
            scale_size_continuous(range = c(2, 15), name = "|Score|") +
            theme_minimal(base_size = 14) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                  axis.text.y = element_text(size = 10),
                  panel.grid.major = element_line(color = "gray90"),
                  panel.border = element_rect(fill = NA, color = "black"),
                  plot.title = element_text(face = "bold", size = 16, hjust = 0.5)) +
            labs(title = "Positive Signature Scores (Enriched Cell Types)",
                 x = "Histopathological Region",
                 y = "Cell Type")
          print(p_pos)
        } else {
          plot.new()
          text(0.5, 0.5, "No positive enrichment detected", adj = 0.5, cex = 1.5)
        }
        
        # Page 5: Negative Enrichment Bubble Plot
        plot_data_neg <- rv$query_results[rv$query_results$Signature_Score < 0, , drop = FALSE]
        
        if (nrow(plot_data_neg) > 0) {
          p_neg <- ggplot(plot_data_neg, aes(x = Region, y = Minor_Cell_Type)) +
            geom_point(aes(size = abs(Signature_Score), color = Signature_Score)) +
            scale_color_gradient2(low = "blue", mid = "white", high = "white",
                                  midpoint = 0, name = "Signature Score") +
            scale_size_continuous(range = c(2, 15), name = "|Score|") +
            theme_minimal(base_size = 14) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                  axis.text.y = element_text(size = 10),
                  panel.grid.major = element_line(color = "gray90"),
                  panel.border = element_rect(fill = NA, color = "black"),
                  plot.title = element_text(face = "bold", size = 16, hjust = 0.5)) +
            labs(title = "Negative Signature Scores (Depleted Cell Types)",
                 x = "Histopathological Region",
                 y = "Cell Type")
          print(p_neg)
        } else {
          plot.new()
          text(0.5, 0.5, "No negative enrichment detected", adj = 0.5, cex = 1.5)
        }
      }
      
      dev.off()
    }
  )
  
  observeEvent(input$reset_analysis, {
    rv$spatial_expr <- NULL
    rv$metadata <- NULL
    rv$scrna_obj <- NULL
    rv$scrna_obj_filtered <- NULL
    rv$scrna_raw_matrix <- NULL
    rv$scrna_raw_metadata <- NULL
    rv$filtered_genes <- NULL
    rv$region_specific_genes <- NULL
    rv$region_gene_list <- NULL
    rv$region_gene_list_overlap <- NULL
    rv$query_results <- NULL
    rv$pca_result <- NULL
    rv$pca_data <- NULL
    rv$var_explained <- NULL
    rv$overlap_matrix <- NULL
    rv$region_gene_list_for_venn <- NULL
    rv$scrna_matrix_loaded <- FALSE
    rv$scrna_genes_loaded <- FALSE
    rv$scrna_barcodes_loaded <- FALSE
    rv$scrna_metadata_loaded <- FALSE
    rv$scrna_matrix_data <- NULL
    rv$scrna_genes_data <- NULL
    rv$scrna_barcodes_data <- NULL
    rv$scrna_metadata_data <- NULL
    rv$proteomics_validated <- FALSE
    rv$scrna_loaded <- FALSE
    rv$scrna_processed <- FALSE
    rv$region_analysis_complete <- FALSE
    rv$query_method_complete <- FALSE
  })
  
  # Download Example Data Handlers
  # Download Example Data Handlers
  output$download_example_expr <- downloadHandler(
    filename = function() {
      "spatial_proteomics_example_data.xlsx"
    },
    content = function(file) {
      # Get the app directory
      app_dir <- "/home/ubuntu/my-hpc-app/data"
      example_file <- file.path(app_dir, "spatial proteomics expression matrix example.xlsx")
      
      if (file.exists(example_file)) {
        file.copy(example_file, file)
      } else {
        stop("Example proteomics data file not found. Please ensure 'spatial proteomics expression matrix example.xlsx' is in the app directory.")
      }
    }
  )
  
  output$download_example_metadata <- downloadHandler(
    filename = function() {
      "spatial_proteomics_metadata_example.xlsx"
    },
    content = function(file) {
      # Get the app directory
      app_dir <- "/home/ubuntu/my-hpc-app/data"
      example_file <- file.path(app_dir, "spatial proteomics metadata example.xlsx")
      
      if (file.exists(example_file)) {
        file.copy(example_file, file)
      } else {
        stop("Example metadata file not found. Please ensure 'spatial proteomics expression matrix example.xlsx' is in the app directory.")
      }
    }
  )
  
  # Output flags for conditional panels
  output$region_analysis_complete <- reactive({ rv$region_analysis_complete })
  output$query_method_complete <- reactive({ rv$query_method_complete })
  
  outputOptions(output, "region_analysis_complete", suspendWhenHidden = FALSE)
  outputOptions(output, "query_method_complete", suspendWhenHidden = FALSE)
  
  # Analysis Summary
  output$analysis_summary <- renderText({
    summary_text <- "Analysis Progress:\n\n"
    
    if (rv$proteomics_validated) {
      summary_text <- paste0(summary_text, "✓ Step 1: Proteomics Data Validated\n")
    }
    if (rv$region_analysis_complete) {
      summary_text <- paste0(summary_text, "✓ Step 2: Region-Specific Proteins Identified\n")
    }
    if (rv$scrna_loaded) {
      summary_text <- paste0(summary_text, "✓ Step 3: scRNA-seq Data Loaded\n")
    }
    if (rv$scrna_processed) {
      summary_text <- paste0(summary_text, "✓ Step 4: scRNA-seq Data Processed\n")
    }
    if (rv$query_method_complete) {
      summary_text <- paste0(summary_text, "✓ Step 5: Query Method Analysis Complete\n")
    }
    
    if (!any(c(rv$proteomics_validated, rv$region_analysis_complete, 
               rv$scrna_loaded, rv$scrna_processed, rv$query_method_complete))) {
      summary_text <- "No analysis completed yet. Please start with Step 1: Proteomics Upload."
    }
    
    summary_text
  })
}

# Run the application
shinyApp(ui = ui, server = server)
