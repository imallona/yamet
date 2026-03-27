## Shared HVF -> scale -> PCA -> UMAP pipeline used by ecker.Rmd and
## crc_embeddings.Rmd. Source after ggtheme.R.
##
## Input convention: mat is features x cells (rows = features, cols = cells).

## Core pipeline: HVF selection, feature scaling, PCA, UMAP.
## Returns list(umap = cells x 2 coords, pca = cells x n_pcs, hvf_idx = row indices)
## or NULL if the data is too small.
run_embedding <- function(mat,
                          n_hvf = 1000L,
                          n_pcs = 50L,
                          n_neighbors = 15L,
                          min_dist = 0.3,
                          seed = 42L) {
  if (is.null(mat) || ncol(mat) < 4L) return(NULL)

  mat[!is.finite(mat)] <- NA

  vars <- matrixStats::rowVars(mat, useNames = FALSE, na.rm = TRUE)
  n_keep <- min(as.integer(n_hvf), sum(!is.na(vars) & vars > 0))
  if (n_keep < 2L) {
    message("run_embedding: too few variable features (", n_keep, "); returning NULL")
    return(NULL)
  }
  hvf_idx <- order(vars, decreasing = TRUE)[seq_len(n_keep)]
  sub_mat <- mat[hvf_idx, , drop = FALSE]

  row_means <- rowMeans(sub_mat, na.rm = TRUE)
  for (i in seq_len(nrow(sub_mat))) {
    nas <- is.na(sub_mat[i, ])
    if (any(nas)) sub_mat[i, nas] <- row_means[i]
  }

  ## scale each feature to mean 0, sd 1 across cells
  ## scale() operates on columns so transpose, scale, transpose back
  scaled <- t(scale(t(sub_mat)))
  bad <- apply(scaled, 1, function(x) all(is.na(x) | is.nan(x)))
  scaled <- scaled[!bad, , drop = FALSE]
  if (nrow(scaled) == 0L) return(NULL)

  n_pcs_use <- min(as.integer(n_pcs), nrow(scaled) - 1L, ncol(scaled) - 1L)
  if (n_pcs_use < 2L) {
    message("run_embedding: too few PCs available (", n_pcs_use, "); returning NULL")
    return(NULL)
  }
  pca_scores <- prcomp(t(scaled), center = FALSE, scale. = FALSE,
                       rank. = n_pcs_use)$x
  message("run_embedding: HVF = ", n_keep, ", PCs = ", n_pcs_use,
          ", cells = ", ncol(mat))

  nn <- min(as.integer(n_neighbors), nrow(pca_scores) - 1L)
  set.seed(seed)
  coords <- uwot::umap(pca_scores, n_neighbors = nn, min_dist = min_dist,
                       n_epochs = 500, verbose = FALSE)

  list(umap = coords, pca = pca_scores, hvf_idx = hvf_idx)
}

## Mean silhouette width for cluster labels on a coordinate matrix.
## Returns NA when fewer than 2 unique labels are present.
sil_score <- function(coords, labels) {
  lv <- as.integer(as.factor(labels))
  if (length(unique(lv)) < 2L) return(NA_real_)
  mean(cluster::silhouette(lv, dist(coords))[, "sil_width"])
}
