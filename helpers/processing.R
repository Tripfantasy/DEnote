# helpers/processing.R
# Seurat pipeline wrappers used throughout the app

library(Seurat)
library(dplyr)

# ── Data Loading ─────────────────────────────────────────────────────────────

load_rds <- function(path) {
  obj <- readRDS(path)
  if (!inherits(obj, "Seurat")) stop("File does not contain a Seurat object.")
  obj
}

load_10x <- function(matrix_path, barcodes_path, features_path,
                     min_cells = 3, min_features = 200) {
  # Write files to a temp dir so Read10X can find them with expected names
  tmp <- tempdir()
  file.copy(matrix_path,   file.path(tmp, "matrix.mtx.gz"),   overwrite = TRUE)
  file.copy(barcodes_path, file.path(tmp, "barcodes.tsv.gz"), overwrite = TRUE)
  file.copy(features_path, file.path(tmp, "features.tsv.gz"), overwrite = TRUE)

  counts <- Read10X(data.dir = tmp)
  # Read10X may return a list for multimodal data; take RNA slot if so
  if (is.list(counts)) counts <- counts[["Gene Expression"]]

  CreateSeuratObject(counts = counts,
                     min.cells = min_cells,
                     min.features = min_features)
}

# ── Standard RNA Pre-processing ──────────────────────────────────────────────

preprocess_seurat <- function(obj, n_pcs = 30, resolution = 0.5,
                              algorithm = 1, progress_fn = NULL) {
  msg <- function(m) if (!is.null(progress_fn)) progress_fn(m)

  msg("Normalizing data...")
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^[Mm][Tt]-")
  obj <- NormalizeData(obj, verbose = FALSE)

  msg("Finding variable features...")
  obj <- FindVariableFeatures(obj, nfeatures = 3000, verbose = FALSE)

  msg("Scaling data...")
  obj <- ScaleData(obj, vars.to.regress = "percent.mt", verbose = FALSE)

  msg("Running PCA...")
  obj <- RunPCA(obj, npcs = n_pcs, verbose = FALSE)

  msg("Running UMAP (2D)...")
  obj <- RunUMAP(obj, dims = seq_len(n_pcs), verbose = FALSE)

  msg("Finding neighbors...")
  obj <- FindNeighbors(obj, dims = seq_len(n_pcs), verbose = FALSE)

  msg("Finding clusters...")
  obj <- FindClusters(obj, resolution = resolution,
                      algorithm = algorithm, verbose = FALSE)
  obj
}

# ── Re-clustering (resolution / algorithm change only) ───────────────────────

recluster <- function(obj, resolution = 0.5, algorithm = 1) {
  # Assumes FindNeighbors has already been run
  FindClusters(obj, resolution = resolution,
               algorithm = algorithm, verbose = FALSE)
}

# ── 3D UMAP ──────────────────────────────────────────────────────────────────

compute_umap3d <- function(obj, n_pcs = 30, reduction = "pca") {
  RunUMAP(obj,
          dims          = seq_len(min(n_pcs, ncol(obj[[reduction]]))),
          n.components  = 3L,
          reduction     = reduction,
          reduction.name = "umap3d",
          verbose       = FALSE)
}

# ── De Novo Marker Analysis ───────────────────────────────────────────────────

find_all_markers <- function(obj, assay = "RNA",
                              only_pos = TRUE,
                              min_pct  = 0.25,
                              logfc    = 0.25,
                              group_by = "seurat_clusters") {
  DefaultAssay(obj) <- assay
  Idents(obj)       <- group_by
  FindAllMarkers(obj,
                 only.pos          = only_pos,
                 min.pct           = min_pct,
                 logfc.threshold   = logfc,
                 verbose           = FALSE)
}

# ── Module Score Labeling (adapted from yihui_cellLabelling.Rmd) ─────────────

run_module_scoring <- function(obj, marker_list, assay = "RNA",
                               threshold = 0.05) {
  DefaultAssay(obj) <- assay

  # Filter out genes not in the dataset
  valid_genes <- rownames(obj)
  marker_list <- lapply(marker_list, function(g) intersect(g, valid_genes))
  marker_list <- Filter(function(g) length(g) > 0, marker_list)

  if (length(marker_list) == 0) stop("None of the provided marker genes were found in the dataset.")

  obj <- AddModuleScore(obj,
                        features = marker_list,
                        name     = names(marker_list),
                        assay    = assay)

  score_cols <- paste0(names(marker_list), seq_along(marker_list))
  obj@meta.data$predicted_celltype <- apply(
    obj@meta.data[, score_cols, drop = FALSE], 1,
    function(x) {
      if (all(is.na(x))) return(NA_character_)
      if (max(x, na.rm = TRUE) < threshold) return("Unknown")
      names(marker_list)[which.max(x)]
    }
  )
  obj
}

# ── Marker list parser ────────────────────────────────────────────────────────
# Parses text like:
#   B-cell: Pax5, Cd79a, Ms4a1
#   Macrophage: F13a1, Cd300a
parse_marker_text <- function(text) {
  lines  <- strsplit(trimws(text), "\n")[[1]]
  lines  <- lines[nzchar(trimws(lines))]
  result <- list()
  for (line in lines) {
    parts <- strsplit(line, ":", fixed = TRUE)[[1]]
    if (length(parts) < 2) next
    cell_type <- trimws(parts[1])
    genes     <- trimws(strsplit(trimws(parts[2]), ",")[[1]])
    genes     <- genes[nzchar(genes)]
    if (length(genes) > 0) result[[cell_type]] <- genes
  }
  result
}

# ── Cluster similarity (Bhattacharyya distance) ───────────────────────────────

compute_cluster_similarity <- function(obj, markers_df,
                                       top_n = 20,
                                       cluster_col = "seurat_clusters") {
  if (is.null(markers_df) || nrow(markers_df) == 0) return(NULL)

  clusters <- as.character(sort(unique(obj@meta.data[[cluster_col]])))

  # Top N markers per cluster
  top_genes <- markers_df |>
    group_by(cluster) |>
    slice_max(order_by = avg_log2FC, n = top_n, with_ties = FALSE) |>
    ungroup()

  all_genes <- unique(top_genes$gene)
  DefaultAssay(obj) <- "RNA"

  # Average expression per cluster
  avg_exp <- AverageExpression(obj,
                               features = all_genes,
                               group.by = cluster_col,
                               layer    = "data",
                               verbose  = FALSE)$RNA
  avg_exp <- as.matrix(avg_exp)

  # Bhattacharyya coefficient between pairs
  n  <- ncol(avg_exp)
  bc <- matrix(0, n, n, dimnames = list(colnames(avg_exp), colnames(avg_exp)))

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      p <- avg_exp[, i]
      q <- avg_exp[, j]
      # Normalise to distributions
      p <- p / (sum(p) + 1e-9)
      q <- q / (sum(q) + 1e-9)
      bc[i, j] <- sum(sqrt(p * q))
    }
  }
  bc
}

# ── Ensure Reductions (compute if missing) ───────────────────────────────────

ensure_pca <- function(obj, n_pcs = 30) {
  if ("pca" %in% names(obj@reductions)) return(obj)
  if (length(VariableFeatures(obj)) == 0)
    obj <- FindVariableFeatures(obj, nfeatures = 3000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = n_pcs, verbose = FALSE)
  obj
}

ensure_umap <- function(obj, n_pcs = 30) {
  obj <- ensure_pca(obj, n_pcs)
  if ("umap" %in% names(obj@reductions)) return(obj)
  avail_dims <- ncol(Embeddings(obj, "pca"))
  obj <- RunUMAP(obj, dims = seq_len(min(n_pcs, avail_dims)), verbose = FALSE)
  obj
}

ensure_tsne <- function(obj, n_pcs = 30) {
  obj <- ensure_pca(obj, n_pcs)
  if ("tsne" %in% names(obj@reductions)) return(obj)
  avail_dims <- ncol(Embeddings(obj, "pca"))
  obj <- RunTSNE(obj, dims = seq_len(min(n_pcs, avail_dims)), verbose = FALSE)
  obj
}

# ── Utility ───────────────────────────────────────────────────────────────────

available_reductions <- function(obj) names(obj@reductions)

active_cluster_col <- function(obj) {
  # Return the most recently created cluster column
  meta_cols  <- colnames(obj@meta.data)
  clust_cols <- grep("_res\\.|seurat_clusters", meta_cols, value = TRUE)
  if (length(clust_cols) > 0) tail(clust_cols, 1) else "seurat_clusters"
}

apply_manual_labels <- function(obj, label_map, cluster_col = "seurat_clusters") {
  clusters <- as.character(obj@meta.data[[cluster_col]])
  obj@meta.data$manual_label <- dplyr::recode(clusters, !!!label_map)
  obj
}
