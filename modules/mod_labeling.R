# modules/mod_labeling.R
library(rlang)

MARKER_PLACEHOLDER <- "B-cell: Pax5, Cd79a, Ms4a1
NK/T-cell: Cd3g, Cd3e, Il2rb
Macrophage: F13a1, Cd300a, Fabp4
Microglia: Siglech, Gpr34, P2ry12
Fibroblast: Col1a1, Col1a2, Dcn
Endothelial: Pecam1, Tie1, Robo4
Pericyte: Rgs5, Rgs4, Notch3
Oligodendrocytes: Mog, Gjb1, Hapln2
Astrocytes: Aldh1l1, Aqp4, Slc1a2
Neuron: Rbfox3, Cacng3, Baiap3"

mod_labeling_ui <- function(id) {
  ns <- NS(id)
  tagList(
    navset_card_underline(
      full_screen = TRUE,

      # Tab 1: Module Score Labeling
      nav_panel(
        "Module Score Labeling",
        layout_columns(
          col_widths = c(4, 8),
          card(
            card_header("Marker Gene Input"),
            p(class = "text-muted small",
              "One cell type per line in the format:",
              code("CellType: Gene1, Gene2, Gene3")),
            textAreaInput(
              ns("marker_text"),
              label       = NULL,
              value       = "",
              placeholder = MARKER_PLACEHOLDER,
              rows        = 14,
              resize      = "vertical"
            ),
            sliderInput(ns("threshold"), "Unknown threshold",
                        0, 0.5, 0.05, step = 0.01),
            selectInput(ns("score_assay"), "Assay for scoring",
                        choices = c("RNA", "SCT"), selected = "RNA"),
            actionButton(ns("score_btn"), "Score & Label",
                         class = "btn-primary w-100", icon = icon("tag"))
          ),
          tagList(
            uiOutput(ns("score_status")),
            navset_card_underline(
              full_screen = TRUE,
              nav_panel(
                tagList(bsicons::bs_icon("grid-3x3-gap"), " Clusters"),
                shinycssloaders::withSpinner(
                  plotlyOutput(ns("cluster_umap"), height = "480px"),
                  type = 6
                )
              ),
              nav_panel(
                tagList(bsicons::bs_icon("tags"), " Predicted Cell Types"),
                shinycssloaders::withSpinner(
                  plotlyOutput(ns("celltype_umap"), height = "480px"),
                  type = 6
                )
              ),
              nav_panel(
                tagList(bsicons::bs_icon("table"), " Score Summary"),
                shinycssloaders::withSpinner(
                  DT::dataTableOutput(ns("score_summary_table")),
                  type = 6
                )
              )
            )
          )
        )
      ),

      # Tab 2: Manual Labeling
      nav_panel(
        "Manual Labeling",
        layout_columns(
          col_widths = c(4, 8),
          card(
            card_header("Assign Labels"),
            p(class = "text-muted small",
              "Edit the label for each cluster below. Pre-populated with module score predictions if available."),
            uiOutput(ns("label_inputs")),
            hr(),
            actionButton(ns("apply_labels_btn"), "Apply Labels",
                         class = "btn-success w-100", icon = icon("check")),
            hr(),
            downloadButton(ns("download_labels"), "Download Labels (.csv)",
                           class = "btn-outline-secondary w-100")
          ),
          card(
            full_screen = TRUE,
            card_header("Labeled UMAP Preview"),
            shinycssloaders::withSpinner(
              plotlyOutput(ns("manual_umap"), height = "500px"),
              type = 6
            )
          )
        )
      )
    )
  )
}

mod_labeling_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Ensure assay selector matches loaded object
    observeEvent(rv$seurat, {
      req(rv$seurat)
      assays <- Assays(rv$seurat)
      updateSelectInput(session, "score_assay",
                        choices  = assays,
                        selected = if ("RNA" %in% assays) "RNA" else assays[1])
    })

    # Module Score Labeling
    observeEvent(input$score_btn, {
      req(rv$seurat, nzchar(input$marker_text))
      tryCatch({
        marker_list <- parse_marker_text(input$marker_text)
        if (length(marker_list) == 0) {
          showNotification("No valid marker entries parsed. Check the format.", type = "warning")
          return()
        }
        withProgress(message = "Running module scoring...", {
          rv$seurat <- run_module_scoring(
            rv$seurat,
            marker_list = marker_list,
            assay       = input$score_assay,
            threshold   = input$threshold
          )
        })
        showNotification("Module scoring complete.", type = "message")
      }, error = function(e) {
        showNotification(paste("Error:", conditionMessage(e)), type = "error")
      })
    })

    output$score_status <- renderUI({
      obj <- rv$seurat
      if (is.null(obj) || !"predicted_celltype" %in% colnames(obj@meta.data)) {
        return(div(class = "alert alert-secondary mb-2",
                   "Provide marker genes and click 'Score & Label' to run module scoring."))
      }
      n_unknown <- sum(obj@meta.data$predicted_celltype == "Unknown", na.rm = TRUE)
      pct       <- round(100 * n_unknown / ncol(obj), 1)
      div(class = "alert alert-success mb-2",
          paste0("Labeling complete. ", n_unknown, " cells (", pct,
                 "%) assigned 'Unknown'."))
    })

    output$cluster_umap <- renderPlotly({
      req(rv$seurat)
      obj <- rv$seurat
      red <- if ("umap" %in% names(obj@reductions)) "umap" else available_reductions(obj)[1]
      umap2d_plotly(obj, color_by = active_cluster_col(obj),
                    reduction = red, title = "Clusters")
    })

    output$celltype_umap <- renderPlotly({
      req(rv$seurat)
      obj <- rv$seurat
      if (!"predicted_celltype" %in% colnames(obj@meta.data)) return(NULL)
      red <- if ("umap" %in% names(obj@reductions)) "umap" else available_reductions(obj)[1]
      umap2d_plotly(obj, color_by = "predicted_celltype",
                    reduction = red, title = "Predicted Cell Types")
    })

    output$score_summary_table <- DT::renderDataTable({
      req(rv$seurat)
      obj <- rv$seurat
      if (!"predicted_celltype" %in% colnames(obj@meta.data)) return(NULL)

      cluster_col <- active_cluster_col(obj)
      df <- obj@meta.data |>
        dplyr::group_by(!!sym(cluster_col), predicted_celltype) |>
        dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") |>
        dplyr::group_by(!!sym(cluster_col)) |>
        dplyr::mutate(pct = round(100 * n_cells / sum(n_cells), 1)) |>
        dplyr::arrange(!!sym(cluster_col), dplyr::desc(n_cells))

      DT::datatable(df, rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE),
                    class = "table table-sm table-striped")
    })

    # Manual Label Inputs
    output$label_inputs <- renderUI({
      req(rv$seurat)
      obj         <- rv$seurat
      cluster_col <- active_cluster_col(obj)
      clusters    <- sort(unique(as.character(obj@meta.data[[cluster_col]])))

      # Pre-populate with module score if available
      has_module <- "predicted_celltype" %in% colnames(obj@meta.data)
      defaults <- if (has_module) {
        obj@meta.data |>
          dplyr::group_by(!!sym(cluster_col), predicted_celltype) |>
          dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
          dplyr::group_by(!!sym(cluster_col)) |>
          dplyr::slice_max(n, n = 1, with_ties = FALSE) |>
          dplyr::select(!!sym(cluster_col), predicted_celltype) |>
          tibble::deframe()
      } else {
        setNames(rep("", length(clusters)), clusters)
      }

      lapply(clusters, function(cl) {
        fluidRow(
          column(3, tags$label(paste0("Cluster ", cl), style = "padding-top:7px;")),
          column(9, textInput(
            ns(paste0("label_", cl)),
            label = NULL,
            value = defaults[[cl]] %||% ""
          ))
        )
      })
    })

    observeEvent(input$apply_labels_btn, {
      req(rv$seurat)
      obj         <- rv$seurat
      cluster_col <- active_cluster_col(obj)
      clusters    <- sort(unique(as.character(obj@meta.data[[cluster_col]])))

      label_map <- sapply(clusters, function(cl) {
        val <- input[[paste0("label_", cl)]]
        if (is.null(val) || !nzchar(val)) cl else val
      })

      rv$seurat <- apply_manual_labels(obj, label_map, cluster_col = cluster_col)
      rv$labels <- label_map
      showNotification("Labels applied.", type = "message")
    })

    output$manual_umap <- renderPlotly({
      req(rv$seurat)
      obj <- rv$seurat
      red <- if ("umap" %in% names(obj@reductions)) "umap" else available_reductions(obj)[1]

      if ("manual_label" %in% colnames(obj@meta.data)) {
        umap2d_plotly(obj, color_by = "manual_label",
                      reduction = red, title = "Manual Labels")
      } else {
        umap2d_plotly(obj, color_by = active_cluster_col(obj),
                      reduction = red, title = "Clusters (labels not yet applied)")
      }
    })

    output$download_labels <- downloadHandler(
      filename = function() paste0("cell_labels_", Sys.Date(), ".csv"),
      content  = function(file) {
        req(rv$seurat)
        obj         <- rv$seurat
        cluster_col <- active_cluster_col(obj)
        df <- data.frame(
          barcode     = colnames(obj),
          cluster     = as.character(obj@meta.data[[cluster_col]]),
          manual_label = if ("manual_label" %in% colnames(obj@meta.data))
                           obj@meta.data$manual_label
                         else NA_character_,
          predicted_celltype = if ("predicted_celltype" %in% colnames(obj@meta.data))
                                  obj@meta.data$predicted_celltype
                                else NA_character_
        )
        write.csv(df, file, row.names = FALSE)
      }
    )

  })
}

# Null operator
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
