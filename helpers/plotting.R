# helpers/plotting.R
# Shared plotting utilities: plotly 2D/3D UMAPs, violin plots, heatmaps

library(ggplot2)
library(plotly)
library(Seurat)
library(dplyr)

# ── Palette ───────────────────────────────────────────────────────────────────

cluster_palette <- function(n) {
  if (n <= 8) {
    RColorBrewer::brewer.pal(max(3, n), "Set2")[seq_len(n)]
  } else if (n <= 12) {
    RColorBrewer::brewer.pal(n, "Set3")
  } else {
    colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n)
  }
}

# ── 2D UMAP (plotly) ─────────────────────────────────────────────────────────

umap2d_plotly <- function(obj, color_by = "seurat_clusters",
                           reduction = "umap", title = NULL,
                           point_size = 3, alpha = 0.7) {
  emb  <- as.data.frame(Embeddings(obj, reduction))
  red_label <- toupper(reduction)
  colnames(emb)[1:2] <- c("Dim_1", "Dim_2")
  emb$cell_id <- rownames(emb)

  # Attach colour variable from metadata
  if (color_by %in% colnames(obj@meta.data)) {
    emb$color_var <- as.character(obj@meta.data[[color_by]])
  } else {
    emb$color_var <- "All cells"
  }

  n_levels <- length(unique(emb$color_var))
  pal      <- cluster_palette(n_levels)
  names(pal) <- sort(unique(emb$color_var))

  p <- plot_ly(
    data   = emb,
    x      = ~Dim_1,
    y      = ~Dim_2,
    color  = ~color_var,
    colors = pal,
    type   = "scatter",
    mode   = "markers",
    text   = ~paste0("Cell: ", cell_id, "<br>", color_by, ": ", color_var),
    hoverinfo = "text",
    marker = list(size = point_size, opacity = alpha)
  ) |>
    layout(
      title  = list(text = title, font = list(size = 14)),
      xaxis  = list(title = paste0(red_label, " 1")),
      yaxis  = list(title = paste0(red_label, " 2")),
      legend = list(title = list(text = color_by)),
      dragmode = "select"
    )
  p
}

# ── 2D UMAP with continuous expression overlay ────────────────────────────────

umap2d_feature_plotly <- function(obj, gene, reduction = "umap",
                                   point_size = 3, alpha = 0.7) {
  emb <- as.data.frame(Embeddings(obj, reduction))
  red_label <- toupper(reduction)
  colnames(emb)[1:2] <- c("Dim_1", "Dim_2")
  emb$cell_id <- rownames(emb)

  DefaultAssay(obj) <- "RNA"
  if (!gene %in% rownames(obj)) {
    stop(paste0("Gene '", gene, "' not found in dataset."))
  }
  expr <- as.numeric(GetAssayData(obj, layer = "data")[gene, ])
  emb$expression <- expr

  plot_ly(
    data  = emb,
    x     = ~Dim_1,
    y     = ~Dim_2,
    color = ~expression,
    colors = c("#d3d3d3", "#fee08b", "#d73027"),
    type  = "scatter",
    mode  = "markers",
    text  = ~paste0("Cell: ", cell_id, "<br>Expression: ", round(expression, 3)),
    hoverinfo = "text",
    marker = list(size = point_size, opacity = alpha)
  ) |>
    layout(
      title = list(text = paste0(gene, " — Expression"), font = list(size = 14)),
      xaxis = list(title = paste0(red_label, " 1")),
      yaxis = list(title = paste0(red_label, " 2"))
    ) |>
    colorbar(title = "Log-norm\nExpression")
}

# ── 3D UMAP (plotly) ─────────────────────────────────────────────────────────

umap3d_plotly <- function(obj, color_by = "seurat_clusters",
                           reduction = "umap3d", title = "3D UMAP",
                           point_size = 2, alpha = 0.7) {
  if (!reduction %in% names(obj@reductions)) return(NULL)

  emb  <- as.data.frame(Embeddings(obj, reduction))
  colnames(emb) <- c("UMAP3D_1", "UMAP3D_2", "UMAP3D_3")
  emb$cell_id   <- rownames(emb)

  if (color_by %in% colnames(obj@meta.data)) {
    emb$color_var <- as.character(obj@meta.data[[color_by]])
  } else {
    emb$color_var <- "All cells"
  }

  n_levels <- length(unique(emb$color_var))
  pal      <- cluster_palette(n_levels)
  names(pal) <- sort(unique(emb$color_var))

  plot_ly(
    data   = emb,
    x      = ~UMAP3D_1,
    y      = ~UMAP3D_2,
    z      = ~UMAP3D_3,
    color  = ~color_var,
    colors = pal,
    type   = "scatter3d",
    mode   = "markers",
    text   = ~paste0("Cell: ", cell_id, "<br>", color_by, ": ", color_var),
    hoverinfo = "text",
    marker = list(size = point_size, opacity = alpha)
  ) |>
    layout(
      title  = list(text = title, font = list(size = 14)),
      scene  = list(
        xaxis = list(title = "UMAP3D 1"),
        yaxis = list(title = "UMAP3D 2"),
        zaxis = list(title = "UMAP3D 3")
      ),
      legend = list(title = list(text = color_by))
    )
}

# ── Violin plot ───────────────────────────────────────────────────────────────

violin_plot <- function(obj, gene, cluster_col = "seurat_clusters") {
  DefaultAssay(obj) <- "RNA"
  if (!gene %in% rownames(obj)) return(NULL)

  expr <- as.numeric(GetAssayData(obj, layer = "data")[gene, ])
  df   <- data.frame(
    expression = expr,
    cluster    = as.character(obj@meta.data[[cluster_col]])
  )
  n_clusters <- length(unique(df$cluster))
  pal        <- cluster_palette(n_clusters)

  ggplot(df, aes(x = cluster, y = expression, fill = cluster)) +
    geom_violin(scale = "width", trim = TRUE, color = NA) +
    geom_jitter(width = 0.2, size = 0.4, alpha = 0.3, color = "grey30") +
    scale_fill_manual(values = pal) +
    labs(title = paste0(gene, " — expression per cluster"),
         x = "Cluster", y = "Log-normalised expression") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none")
}

# ── Cluster similarity heatmap ────────────────────────────────────────────────

similarity_heatmap <- function(bc_matrix) {
  if (is.null(bc_matrix)) return(NULL)

  df <- as.data.frame(bc_matrix) |>
    tibble::rownames_to_column("Cluster_A") |>
    tidyr::pivot_longer(-Cluster_A, names_to = "Cluster_B", values_to = "BC")

  ggplot(df, aes(x = Cluster_B, y = Cluster_A, fill = BC)) +
    geom_tile(color = "white") +
    scale_fill_distiller(
      palette  = "RdYlBu",
      direction = -1,
      limits   = c(0, 1),
      name     = "Bhattacharyya\nCoefficient"
    ) +
    labs(
      title   = "Cluster Similarity (marker expression)",
      subtitle = "High values → similar marker profiles → candidate merge",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ── Static labelled UMAP (ggplot2, for reports) ───────────────────────────────

umap_ggplot <- function(obj, color_by = "seurat_clusters",
                         label = TRUE, reduction = "umap",
                         title = NULL) {
  emb <- as.data.frame(Embeddings(obj, reduction))
  red_label <- toupper(reduction)
  colnames(emb)[1:2] <- c("Dim_1", "Dim_2")
  emb$color_var <- as.character(obj@meta.data[[color_by]])

  centroids <- emb |>
    group_by(color_var) |>
    summarise(Dim_1 = mean(Dim_1), Dim_2 = mean(Dim_2))

  n_levels <- length(unique(emb$color_var))
  pal      <- cluster_palette(n_levels)

  p <- ggplot(emb, aes(x = Dim_1, y = Dim_2, color = color_var)) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_color_manual(values = pal, name = color_by) +
    labs(title = title,
         x = paste0(red_label, " 1"),
         y = paste0(red_label, " 2")) +
    theme_minimal(base_size = 12)

  if (label) {
    p <- p + ggrepel::geom_label_repel(
      data        = centroids,
      aes(label = color_var),
      color       = "black",
      fill        = "white",
      size        = 3,
      label.size  = 0.2,
      label.padding = unit(0.15, "lines"),
      show.legend = FALSE
    )
  }
  p
}
