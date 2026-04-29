# modules/mod_report.R

mod_report_ui <- function(id) {
  ns <- NS(id)
  tagList(
    layout_columns(
      col_widths = c(4, 8),

      # Export controls
      card(
        card_header("Export Options"),

        h6("HTML Report"),
        uiOutput(ns("report_status")),
        actionButton(ns("generate_btn"), "Generate Report",
                     class = "btn-primary w-100 mb-2",
                     icon = icon("file-code")),
        uiOutput(ns("download_ui")),

        hr(),
        h6("Data Tables"),
        downloadButton(ns("download_markers"), "Download Marker Table (.csv)",
                       class = "btn-outline-secondary w-100 mb-2"),
        downloadButton(ns("download_metadata"), "Download Cell Metadata (.csv)",
                       class = "btn-outline-secondary w-100")
      ),

      # Live preview
      card(
        full_screen = TRUE,
        card_header("Report Preview"),
        uiOutput(ns("preview_header")),
        plotOutput(ns("preview_umap"), height = "350px"),
        uiOutput(ns("preview_table_header")),
        DT::dataTableOutput(ns("preview_table"))
      )
    )
  )
}

mod_report_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {

    # Track analysis version for stale detection
    analysis_version <- reactiveVal(0L)
    observeEvent(list(rv$seurat, rv$markers, rv$labels), {
      analysis_version(analysis_version() + 1L)
    }, ignoreNULL = FALSE)

    report_path    <- reactiveVal(NULL)   # path to last rendered HTML
    report_version <- reactiveVal(NA_integer_)  # analysis_version at render time

    report_stale <- reactive({
      is.null(report_path()) ||
        is.na(report_version()) ||
        report_version() < analysis_version()
    })

    # Report status badge
    output$report_status <- renderUI({
      if (is.null(rv$seurat)) {
        div(class = "alert alert-secondary small py-1 mb-2",
            "Load data to enable report generation.")
      } else if (is.null(report_path())) {
        div(class = "alert alert-secondary small py-1 mb-2",
            "Report not yet generated.")
      } else if (report_stale()) {
        div(class = "alert alert-warning small py-1 mb-2",
            bsicons::bs_icon("exclamation-triangle-fill"), " ",
            "Analysis has changed — regenerate to update.")
      } else {
        div(class = "alert alert-success small py-1 mb-2",
            bsicons::bs_icon("check-circle-fill"), " ",
            "Report is up to date.")
      }
    })

    # Show download button only when a report exists
    output$download_ui <- renderUI({
      req(report_path())
      downloadButton(session$ns("download_report"), "Download HTML Report",
                     class = "btn-outline-primary w-100 mb-2")
    })

    # Generate report (pre-render to temp file)
    observeEvent(input$generate_btn, {
      req(rv$seurat)

      ver <- analysis_version()
      setProgress(0.1, detail = "Preparing data")

      tmp_html <- tempfile(fileext = ".html")

      withProgress(message = "Rendering HTML report...", value = 0, {
        setProgress(0.1, detail = "Preparing data")

        tmp_html <- tempfile(fileext = ".html")

        # Inject objects directly 
        render_env           <- new.env(parent = globalenv())
        render_env$obj       <- rv$seurat
        render_env$markers   <- rv$markers

        # Locate template relative to app directory
        template <- file.path(
          dirname(rstudioapi::getSourceEditorContext()$path %||%
                    getwd()),
          "report_template.Rmd"
        )
        if (!file.exists(template)) template <- "report_template.Rmd"

        setProgress(0.2, detail = "Running knitr")

        tryCatch({
          rmarkdown::render(
            input       = template,
            output_file = tmp_html,
            params      = list(labels = rv$labels),
            envir       = render_env,
            quiet       = TRUE
          )
          report_path(tmp_html)
          report_version(ver)
          showNotification("Report ready for download.", type = "message")
        }, error = function(e) {
          showNotification(paste("Render failed:", conditionMessage(e)),
                           type = "error", duration = 10)
        })
      })
    })

    # Download: just copy the pre-rendered file
    output$download_report <- downloadHandler(
      filename = function() paste0("cell_labelling_report_", Sys.Date(), ".html"),
      content  = function(file) {
        req(report_path())
        file.copy(report_path(), file, overwrite = TRUE)
      }
    )

    # Live preview data (auto-updates with analysis state)
    preview_data <- reactive({
      obj <- rv$seurat
      if (is.null(obj)) return(NULL)

      cluster_col <- active_cluster_col(obj)
      has_manual  <- "manual_label"       %in% colnames(obj@meta.data)
      has_module  <- "predicted_celltype" %in% colnames(obj@meta.data)
      has_markers <- !is.null(rv$markers) && nrow(rv$markers) > 0

      tbl             <- table(obj@meta.data[[cluster_col]])
      cluster_summary <- data.frame(
        Cluster = names(tbl),
        N_Cells = as.integer(tbl),
        stringsAsFactors = FALSE
      )

      if (has_manual) {
        uniq <- tapply(obj@meta.data$manual_label,
                       obj@meta.data[[cluster_col]],
                       function(x) names(sort(table(x), decreasing = TRUE))[1])
        cluster_summary$Manual_Label <- uniq[cluster_summary$Cluster]
      }
      if (has_module) {
        uniq <- tapply(obj@meta.data$predicted_celltype,
                       obj@meta.data[[cluster_col]],
                       function(x) names(sort(table(x), decreasing = TRUE))[1])
        cluster_summary$Module_Score_Label <- uniq[cluster_summary$Cluster]
      }
      if (has_markers) {
        top5 <- rv$markers |>
          dplyr::group_by(cluster) |>
          dplyr::slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) |>
          dplyr::summarise(Top_Markers = paste(gene, collapse = ", "),
                           .groups = "drop")
        cluster_summary <- dplyr::left_join(
          cluster_summary, top5, by = c("Cluster" = "cluster")
        )
      }

      list(
        obj             = obj,
        cluster_col     = cluster_col,
        has_manual      = has_manual,
        n_cells         = ncol(obj),
        n_clusters      = length(tbl),
        cluster_summary = cluster_summary
      )
    })

    # Preview: header value boxes
    output$preview_header <- renderUI({
      pd <- preview_data()
      if (is.null(pd)) {
        return(div(class = "alert alert-secondary",
                   "Load data to see a live summary here."))
      }
      tagList(
        layout_column_wrap(
          width = "150px", fill = FALSE,
          value_box("Cells",    scales::comma(pd$n_cells), theme = "primary"),
          value_box("Clusters", pd$n_clusters,             theme = "secondary"),
          value_box("Labeled",
                    if (pd$has_manual) "Yes" else "No",
                    theme = if (pd$has_manual) "success" else "warning")
        ),
        tags$h6("Labeled UMAP", class = "mt-3 mb-1")
      )
    })

    output$preview_umap <- renderPlot({
      pd <- preview_data()
      if (is.null(pd)) return(NULL)
      obj      <- pd$obj
      color_by <- if (pd$has_manual) "manual_label" else pd$cluster_col
      red      <- if ("umap" %in% names(obj@reductions)) "umap" else available_reductions(obj)[1]
      umap_ggplot(obj, color_by = color_by, reduction = red,
                  title = "Labeled Clusters")
    })

    output$preview_table_header <- renderUI({
      if (is.null(preview_data())) return(NULL)
      tags$h6("Cluster Summary", class = "mt-3 mb-1")
    })

    output$preview_table <- DT::renderDataTable({
      pd <- preview_data()
      if (is.null(pd)) return(NULL)
      df <- pd$cluster_summary
      colnames(df) <- gsub("_", " ", colnames(df))
      DT::datatable(df, rownames = FALSE,
                    options = list(pageLength = 20, dom = "t"),
                    class = "table table-sm table-striped")
    })

    # Marker table download
    output$download_markers <- downloadHandler(
      filename = function() paste0("markers_", Sys.Date(), ".csv"),
      content  = function(file) {
        req(rv$markers)
        write.csv(rv$markers, file, row.names = FALSE)
      }
    )

    # Cell metadata download
    output$download_metadata <- downloadHandler(
      filename = function() paste0("cell_metadata_", Sys.Date(), ".csv"),
      content  = function(file) {
        req(rv$seurat)
        write.csv(rv$seurat@meta.data, file, row.names = TRUE)
      }
    )

  })
}
