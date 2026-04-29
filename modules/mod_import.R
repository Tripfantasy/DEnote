# modules/mod_import.R

mod_import_ui <- function(id) {
  ns <- NS(id)
  tagList(
    layout_columns(
      col_widths = c(5, 7),
      # Input Panel
      card(
        card_header("Load Data"),
        radioButtons(
          ns("input_type"), "Input format",
          choices  = c("Seurat RDS" = "rds", "10x Feature/Barcode Matrix" = "tenx"),
          selected = "rds"
        ),
        # RDS input
        conditionalPanel(
          condition = sprintf("input['%s'] == 'rds'", ns("input_type")),
          fileInput(ns("rds_file"), "Upload .rds file",
                    accept = ".rds", placeholder = "No file selected")
        ),
        # 10x input
        conditionalPanel(
          condition = sprintf("input['%s'] == 'tenx'", ns("input_type")),
          fileInput(ns("matrix_file"),   "matrix.mtx or matrix.mtx.gz",
                    accept = c(".mtx", ".gz")),
          fileInput(ns("barcodes_file"), "barcodes.tsv or barcodes.tsv.gz",
                    accept = c(".tsv", ".gz")),
          fileInput(ns("features_file"), "features.tsv or features.tsv.gz",
                    accept = c(".tsv", ".gz")),
          accordion(
            open = FALSE,
            accordion_panel(
              "Pre-processing options",
              sliderInput(ns("min_cells"),    "Min cells per gene",    1, 20,  3, step = 1),
              sliderInput(ns("min_features"), "Min features per cell", 50, 500, 200, step = 50),
              sliderInput(ns("n_pcs"),        "Number of PCs",         10, 50, 30, step = 5),
              sliderInput(ns("init_res"),     "Initial resolution",    0.1, 2, 0.5, step = 0.1),
              selectInput(ns("algorithm"),    "Clustering algorithm",
                          choices = c("Louvain" = 1, "Leiden" = 3), selected = 1)
            )
          ),
          actionButton(ns("process_btn"), "Process Data",
                       class = "btn-primary", icon = icon("play"))
        )
      ),

      # QC Summary
      card(
        card_header("Dataset Summary"),
        uiOutput(ns("qc_summary")),
        uiOutput(ns("status_msg"))
      )
    )
  )
}

mod_import_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {

    # Status messages
    status <- reactiveVal(NULL)

    output$status_msg <- renderUI({
      msg <- status()
      if (is.null(msg)) return(NULL)
      class <- if (grepl("^Error", msg)) "alert alert-danger" else "alert alert-info"
      div(class = class, style = "margin-top:10px;", msg)
    })

    # Load RDS
    observeEvent(input$rds_file, {
      req(input$rds_file)
      status("Loading Seurat object...")
      tryCatch({
        obj        <- load_rds(input$rds_file$datapath)
        rv$seurat  <- obj
        status(NULL)
      }, error = function(e) status(paste("Error:", conditionMessage(e))))
    })

    # Process 10x
    observeEvent(input$process_btn, {
      req(input$matrix_file, input$barcodes_file, input$features_file)
      status("Loading count matrix...")
      tryCatch({
        obj <- load_10x(
          matrix_path   = input$matrix_file$datapath,
          barcodes_path = input$barcodes_file$datapath,
          features_path = input$features_file$datapath,
          min_cells     = input$min_cells,
          min_features  = input$min_features
        )
        obj <- withProgress(message = "Processing...", value = 0, {
          preprocess_seurat(
            obj,
            n_pcs      = input$n_pcs,
            resolution = input$init_res,
            algorithm  = as.integer(input$algorithm),
            progress_fn = function(m) {
              incProgress(1/7, detail = m)
            }
          )
        })
        rv$seurat <- obj
        status(NULL)
      }, error = function(e) status(paste("Error:", conditionMessage(e))))
    })

    # QC Summary
    output$qc_summary <- renderUI({
      obj <- rv$seurat
      if (is.null(obj)) {
        return(p(class = "text-muted", "Load a dataset to see summary information."))
      }

      cluster_col <- active_cluster_col(obj)
      n_clusters  <- length(unique(obj@meta.data[[cluster_col]]))
      assays      <- paste(Assays(obj), collapse = ", ")
      reductions  <- paste(names(obj@reductions), collapse = ", ")
      has_3d      <- "umap3d" %in% names(obj@reductions)

      tagList(
        layout_column_wrap(
          width = "120px", fill = FALSE,
          value_box("Cells",        scales::comma(ncol(obj)),    theme = "primary"),
          value_box("Genes",        scales::comma(nrow(obj)),    theme = "secondary"),
          value_box("Clusters",     n_clusters,                  theme = "info"),
          value_box("3D UMAP",      if (has_3d) "Yes" else "No", theme = if (has_3d) "success" else "light")
        ),
        tags$table(
          class = "table table-sm mt-3",
          tags$tr(tags$th("Assays"),     tags$td(assays)),
          tags$tr(tags$th("Reductions"), tags$td(reductions))
        )
      )
    })

  })
}
