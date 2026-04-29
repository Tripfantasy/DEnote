# modules/mod_markers.R

mod_markers_ui <- function(id) {
  ns <- NS(id)
  tagList(
    layout_columns(
      col_widths = c(3, 9),
      # Sidebar
      card(
        card_header("Marker Controls"),
        h6("Feature Expression"),
        textInput(ns("gene_input"), "Gene name", placeholder = "e.g. Pax5"),
        actionButton(ns("plot_gene_btn"), "Plot expression",
                     class = "btn-primary w-100", icon = icon("chart-area")),
        hr(),
        h6("De Novo Marker Analysis"),
        sliderInput(ns("min_pct"),  "Min pct expressing",  0.05, 0.9, 0.25, step = 0.05),
        sliderInput(ns("logfc"),    "Min log2 FC",          0.1, 2.0, 0.25, step = 0.05),
        actionButton(ns("find_markers_btn"), "Find All Markers",
                     class = "btn-warning w-100", icon = icon("magnifying-glass")),
        hr(),
        h6("Cluster Boundary Optimisation"),
        sliderInput(ns("top_n_genes"), "Top N markers for similarity", 5, 50, 20, step = 5),
        sliderInput(ns("bc_threshold"), "Merge threshold (BC)", 0.50, 0.99, 0.90, step = 0.01),
        actionButton(ns("similarity_btn"), "Compute Similarity",
                     class = "btn-outline-secondary w-100",
                     icon = icon("diagram-project"))
      ),

      # Main panel
      navset_card_underline(
        full_screen = TRUE,
        nav_panel(
          "Feature Plot",
          shinycssloaders::withSpinner(
            plotlyOutput(ns("feature_plot"), height = "450px"),
            type = 6, color = "#0d6efd"
          )
        ),
        nav_panel(
          "Violin Plot",
          selectizeInput(ns("violin_gene"), "Gene for violin",
                         choices = NULL, width = "300px",
                         options = list(maxOptions = 10000)),
          shinycssloaders::withSpinner(
            plotOutput(ns("violin_plot"), height = "400px"),
            type = 6, color = "#0d6efd"
          )
        ),
        nav_panel(
          "De Novo Markers",
          uiOutput(ns("markers_status")),
          shinycssloaders::withSpinner(
            DT::dataTableOutput(ns("markers_table")),
            type = 6, color = "#0d6efd"
          )
        ),
        nav_panel(
          "Cluster Similarity",
          uiOutput(ns("similarity_status")),
          shinycssloaders::withSpinner(
            plotOutput(ns("similarity_heatmap"), height = "450px"),
            type = 6, color = "#0d6efd"
          ),
          uiOutput(ns("merge_suggestions"))
        )
      )
    )
  )
}

mod_markers_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {

    # Feature plot
    feature_data <- eventReactive(input$plot_gene_btn, {
      req(rv$seurat, nzchar(input$gene_input))
      input$gene_input
    })

    output$feature_plot <- renderPlotly({
      gene <- feature_data()
      req(gene)
      obj  <- rv$seurat
      reduction <- if ("umap" %in% names(obj@reductions)) "umap" else available_reductions(obj)[1]
      tryCatch(
        umap2d_feature_plotly(obj, gene, reduction = reduction),
        error = function(e) {
          showNotification(conditionMessage(e), type = "warning")
          NULL
        }
      )
    })

    # Violin plot gene selector 
    observeEvent(rv$seurat, {
      req(rv$seurat)
      genes <- rownames(rv$seurat)
      updateSelectizeInput(session, "violin_gene", choices = genes,
                           selected = genes[1], server = TRUE)
    })

    observeEvent(rv$markers, {
      req(rv$markers)
      top_genes <- unique(rv$markers$gene)
      updateSelectizeInput(session, "violin_gene", choices = top_genes,
                           selected = top_genes[1], server = TRUE)
    })

    output$violin_plot <- renderPlot({
      req(rv$seurat, input$violin_gene)
      cluster_col <- active_cluster_col(rv$seurat)
      violin_plot(rv$seurat, input$violin_gene, cluster_col = cluster_col)
    })

    # Find All Markers
    observeEvent(input$find_markers_btn, {
      req(rv$seurat)
      withProgress(message = "Finding markers (this may take several minutes)...", {
        tryCatch({
          cluster_col   <- active_cluster_col(rv$seurat)
          rv$markers    <- find_all_markers(
            rv$seurat,
            min_pct  = input$min_pct,
            logfc    = input$logfc,
            group_by = cluster_col
          )
          showNotification(
            paste0("Found ", nrow(rv$markers), " markers across ",
                   length(unique(rv$markers$cluster)), " clusters."),
            type = "message"
          )
        }, error = function(e) {
          showNotification(paste("Failed:", conditionMessage(e)), type = "error")
        })
      })
    })

    output$markers_status <- renderUI({
      if (is.null(rv$markers)) {
        div(class = "alert alert-secondary",
            "Click 'Find All Markers' to run de novo differential expression analysis.")
      } else {
        div(class = "alert alert-success",
            paste0(nrow(rv$markers), " markers found across ",
                   length(unique(rv$markers$cluster)), " clusters."))
      }
    })

    output$markers_table <- DT::renderDataTable({
      req(rv$markers)
      df <- rv$markers |>
        dplyr::select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj) |>
        dplyr::mutate(
          avg_log2FC = round(avg_log2FC, 3),
          pct.1      = round(pct.1, 3),
          pct.2      = round(pct.2, 3),
          p_val_adj  = formatC(p_val_adj, format = "e", digits = 2)
        )
      DT::datatable(
        df,
        filter    = "top",
        rownames  = FALSE,
        options   = list(pageLength = 15, scrollX = TRUE),
        class     = "table table-striped table-hover table-sm"
      )
    })

    # Cluster similarity
    bc_matrix <- eventReactive(input$similarity_btn, {
      req(rv$seurat, rv$markers)
      withProgress(message = "Computing cluster similarity...", {
        tryCatch(
          compute_cluster_similarity(rv$seurat, rv$markers,
                                     top_n = input$top_n_genes,
                                     cluster_col = active_cluster_col(rv$seurat)),
          error = function(e) {
            showNotification(conditionMessage(e), type = "error")
            NULL
          }
        )
      })
    })

    output$similarity_status <- renderUI({
      mat <- bc_matrix()
      if (!isTruthy(mat)) {
        div(class = "alert alert-secondary",
            "Run 'Find All Markers' first, then click 'Compute Similarity'.")
      } else {
        n <- nrow(mat)
        div(class = "alert alert-success",
            bsicons::bs_icon("check-circle-fill"), " ",
            paste0("Similarity computed for ", n, " cluster",
                   if (n != 1) "s", "."))
      }
    })

    output$similarity_heatmap <- renderPlot({
      req(bc_matrix())
      similarity_heatmap(bc_matrix())
    })

    output$merge_suggestions <- renderUI({
      mat <- bc_matrix()
      if (is.null(mat)) return(NULL)

      thresh <- input$bc_threshold
      pairs <- which(mat > thresh & mat < 1, arr.ind = TRUE)
      pairs <- pairs[pairs[, 1] < pairs[, 2], , drop = FALSE]

      if (nrow(pairs) == 0) {
        return(div(class = "alert alert-success mt-2",
                   "No obvious merge candidates detected (all BC < 0.85)."))
      }

      items <- apply(pairs, 1, function(r) {
        a <- rownames(mat)[r[1]]
        b <- colnames(mat)[r[2]]
        bc_val <- round(mat[r[1], r[2]], 3)
        tags$li(paste0("Clusters ", a, " & ", b, " — BC = ", bc_val))
      })

      tagList(
        div(class = "alert alert-warning mt-2",
            strong(paste0("Merge Candidates (BC > ", input$bc_threshold, "):")),
            tags$ul(items)
        )
      )
    })

  })
}
