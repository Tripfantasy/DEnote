# modules/mod_clustering.R

mod_clustering_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Controls sidebar row
    layout_columns(
      col_widths = c(3, 9),
      # Sidebar card
      card(
        card_header("Clustering Controls"),
        sliderInput(ns("resolution"), "Resolution", 0.1, 3.0, 0.5, step = 0.1),
        selectInput(ns("algorithm"), "Algorithm",
                    choices = c("Louvain" = "1", "Leiden" = "3"), selected = "1"),
        selectInput(ns("reduction"), "Reduction", choices = NULL),
        actionButton(ns("recluster_btn"), "Re-cluster",
                     class = "btn-primary w-100", icon = icon("rotate")),
        hr(),
        input_switch(ns("show_3d"), "Show 3D UMAP", value = FALSE),
        conditionalPanel(
          condition = sprintf("input['%s']", ns("show_3d")),
          actionButton(ns("compute_3d_btn"), "Compute 3D UMAP",
                       class = "btn-outline-secondary w-100 mt-1",
                       icon = icon("cube"))
        ),
        hr(),
        selectInput(ns("color_by"), "Colour by",
                    choices = c("seurat_clusters"), selected = "seurat_clusters"),
        hr(),
        accordion(
          open = FALSE,
          accordion_panel(
            title = tagList(bsicons::bs_icon("cpu"), " Compute Reductions"),
            value = "reductions",
            numericInput(ns("n_pcs_redu"), "PCs to use", value = 30, min = 5, max = 100),
            uiOutput(ns("reduction_status")),
            div(
              class = "d-grid gap-1 mt-2",
              actionButton(ns("run_pca_btn"),  "Run PCA",
                           class = "btn-sm btn-outline-secondary w-100",
                           icon = icon("chart-bar")),
              actionButton(ns("run_umap_btn"), "Run UMAP",
                           class = "btn-sm btn-outline-secondary w-100",
                           icon = icon("circle-dot")),
              actionButton(ns("run_tsne_btn"), "Run tSNE",
                           class = "btn-sm btn-outline-secondary w-100",
                           icon = icon("diagram-project"))
            )
          )
        )
      ),

      # Main plots
      tagList(
        layout_column_wrap(
          width = 1,
          fill = FALSE,
          layout_columns(
            col_widths = c(4, 4, 4),
            value_box("Clusters",        uiOutput(ns("n_clusters")),  theme = "primary",   fill = FALSE),
            value_box("Cells",           uiOutput(ns("n_cells")),     theme = "secondary", fill = FALSE),
            value_box("Median / cluster",uiOutput(ns("med_cells")),   theme = "info",      fill = FALSE)
          )
        ),
        navset_card_underline(
          full_screen = TRUE,
          nav_panel("2D UMAP",
            shinycssloaders::withSpinner(
              plotlyOutput(ns("umap2d"), height = "500px"),
              type = 6, color = "#0d6efd"
            )
          ),
          nav_panel("3D UMAP",
            uiOutput(ns("umap3d_ui"))
          )
        )
      )
    )
  )
}

mod_clustering_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

  
    # Reduction selecter updates
    refresh_reduction_selector <- function(obj) {
      reductions <- grep("umap|pca|tsne|lsi|svd", available_reductions(obj),
                         ignore.case = TRUE, value = TRUE)
      reductions <- reductions[!grepl("3d", reductions, ignore.case = TRUE)]
      selected <- if ("umap" %in% reductions) "umap" else reductions[1]
      updateSelectInput(session, "reduction", choices = reductions, selected = selected)
    }

    observeEvent(rv$seurat, {
      req(rv$seurat)
      obj       <- rv$seurat
      refresh_reduction_selector(obj)

      meta_cols <- colnames(obj@meta.data)
      color_options <- c(
        meta_cols[grep("cluster|ident|label|celltype|predicted",
                       meta_cols, ignore.case = TRUE)]
      )
      color_options <- unique(c("seurat_clusters", color_options))
      color_options <- intersect(color_options, meta_cols)
      updateSelectInput(session, "color_by",
                        choices  = color_options,
                        selected = active_cluster_col(obj))
    })

    # Re-cluster
    observeEvent(input$recluster_btn, {
      req(rv$seurat)
      withProgress(message = "Re-clustering...", {
        tryCatch({
          rv$seurat <- recluster(rv$seurat,
                                 resolution = input$resolution,
                                 algorithm  = as.integer(input$algorithm))
          # Refresh color palette 
          obj <- rv$seurat
          meta_cols     <- colnames(obj@meta.data)
          color_options <- meta_cols[grep("cluster|ident|label|celltype|predicted",
                                          meta_cols, ignore.case = TRUE)]
          color_options <- unique(c("seurat_clusters", color_options))
          color_options <- intersect(color_options, meta_cols)
          updateSelectInput(session, "color_by",
                            choices  = color_options,
                            selected = active_cluster_col(obj))
        }, error = function(e) {
          showNotification(paste("Re-cluster failed:", conditionMessage(e)),
                           type = "error")
        })
      })
    })

    # Compute 3D UMAP
    observeEvent(input$compute_3d_btn, {
      req(rv$seurat)
      withProgress(message = "Computing 3D UMAP (may take a moment)...", {
        tryCatch({
          rv$seurat <- compute_umap3d(rv$seurat)
          showNotification("3D UMAP computed.", type = "message")
        }, error = function(e) {
          showNotification(paste("3D UMAP failed:", conditionMessage(e)), type = "error")
        })
      })
    })

    # Reduction status display
    output$reduction_status <- renderUI({
      req(rv$seurat)
      reds <- names(rv$seurat@reductions)
      targets <- c("pca" = "PCA", "umap" = "UMAP", "tsne" = "tSNE")
      items <- lapply(names(targets), function(r) {
        present <- r %in% reds
        tags$div(
          class = paste0("d-flex align-items-center gap-1 small ",
                         if (present) "text-success" else "text-secondary"),
          bsicons::bs_icon(if (present) "check-circle-fill" else "circle"),
          targets[[r]]
        )
      })
      div(class = "mb-1", do.call(tagList, items))
    })

    # ── Compute PCA ──────────────────────────────────────────────────────
    observeEvent(input$run_pca_btn, {
      req(rv$seurat)
      withProgress(message = "Computing PCA...", {
        tryCatch({
          rv$seurat <- ensure_pca(rv$seurat, n_pcs = input$n_pcs_redu)
          refresh_reduction_selector(rv$seurat)
          showNotification("PCA ready.", type = "message")
        }, error = function(e) {
          showNotification(paste("PCA failed:", conditionMessage(e)), type = "error")
        })
      })
    })

    # Compute UMAP
    observeEvent(input$run_umap_btn, {
      req(rv$seurat)
      withProgress(message = "Computing UMAP...", {
        tryCatch({
          rv$seurat <- ensure_umap(rv$seurat, n_pcs = input$n_pcs_redu)
          refresh_reduction_selector(rv$seurat)
          showNotification("UMAP ready.", type = "message")
        }, error = function(e) {
          showNotification(paste("UMAP failed:", conditionMessage(e)), type = "error")
        })
      })
    })

    # Compute tSNE
    observeEvent(input$run_tsne_btn, {
      req(rv$seurat)
      withProgress(message = "Computing tSNE (may take a moment)...", {
        tryCatch({
          rv$seurat <- ensure_tsne(rv$seurat, n_pcs = input$n_pcs_redu)
          refresh_reduction_selector(rv$seurat)
          showNotification("tSNE ready.", type = "message")
        }, error = function(e) {
          showNotification(paste("tSNE failed:", conditionMessage(e)), type = "error")
        })
      })
    })

    # Value boxes
    cluster_stats <- reactive({
      req(rv$seurat)
      obj         <- rv$seurat
      cluster_col <- active_cluster_col(obj)
      tbl         <- table(obj@meta.data[[cluster_col]])
      list(
        n_clusters = length(tbl),
        n_cells    = ncol(obj),
        med_cells  = median(as.integer(tbl))
      )
    })
    output$n_clusters <- renderUI(strong(cluster_stats()$n_clusters))
    output$n_cells    <- renderUI(strong(scales::comma(cluster_stats()$n_cells)))
    output$med_cells  <- renderUI(strong(scales::comma(cluster_stats()$med_cells)))

    # 2D UMAP
    output$umap2d <- renderPlotly({
      req(rv$seurat, input$reduction, input$color_by)
      obj  <- rv$seurat
      red  <- input$reduction
      cb   <- input$color_by
      if (!red %in% names(obj@reductions)) return(NULL)
      if (!cb  %in% colnames(obj@meta.data)) cb <- active_cluster_col(obj)
      umap2d_plotly(obj, color_by = cb, reduction = red,
                    title = paste0("2D ", toupper(red)))
    })

    # 3D UMAP
    output$umap3d_ui <- renderUI({
      if (!isTRUE(input$show_3d)) {
        return(p(class = "text-muted p-3",
                 "Enable '3D UMAP' toggle and click 'Compute 3D UMAP' to generate."))
      }
      if (!"umap3d" %in% names(rv$seurat@reductions)) {
        return(p(class = "text-muted p-3",
                 "Click 'Compute 3D UMAP' to generate the embedding."))
      }
      shinycssloaders::withSpinner(
        plotlyOutput(ns("umap3d_plot"), height = "500px"),
        type = 6, color = "#0d6efd"
      )
    })

    output$umap3d_plot <- renderPlotly({
      req(rv$seurat, input$color_by)
      obj <- rv$seurat
      cb  <- input$color_by
      if (!cb %in% colnames(obj@meta.data)) cb <- active_cluster_col(obj)
      umap3d_plotly(obj, color_by = cb, title = "3D UMAP")
    })

    # Expose selected cells from 2D UMAP click events
    observeEvent(event_data("plotly_selected", source = "umap2d"), {
      sel <- event_data("plotly_selected", source = "umap2d")
      if (!is.null(sel) && nrow(sel) > 0) {
        rv$selected_cells <- sel$text  # hover text contains cell IDs
      }
    })

  })
}
