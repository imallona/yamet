## Shared embedding utilities used by ecker.Rmd, argelaguet.Rmd, and
## crc_embeddings.Rmd. Source after ggtheme.R.
##
## Convention for run_embedding / run_pca_mat: mat is features x cells.
## Convention for run_umap_wide / regress_out_meth: cells x features (wide df).
## Requires dplyr, uwot, matrixStats, cluster.

## HVF selection, feature scaling, PCA, UMAP.
## Returns list(umap = cells x 2, pca = cells x n_pcs, hvf_idx = row indices)
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

## PCA on a features x cells matrix. Returns cells x PCs score matrix or NULL.
run_pca_mat <- function(mat, n_pcs = 10L) {
  mat[!is.finite(mat)] <- NA
  row_means <- rowMeans(mat, na.rm = TRUE)
  for (i in seq_len(nrow(mat))) {
    nas <- is.na(mat[i, ])
    if (any(nas)) mat[i, nas] <- row_means[i]
  }
  row_sd <- apply(mat, 1, sd, na.rm = TRUE)
  bad <- is.na(row_sd) | row_sd == 0
  mat <- mat[!bad, , drop = FALSE]
  if (nrow(mat) < 2L || ncol(mat) < 2L) return(NULL)
  n_use <- min(as.integer(n_pcs), nrow(mat) - 1L, ncol(mat) - 1L)
  if (n_use < 1L) return(NULL)
  prcomp(t(mat), center = TRUE, scale. = TRUE, rank. = n_use)$x
}

## UMAP on a compact cells x dims score matrix (no HVF/PCA needed).
run_umap_scores <- function(scores, seed = 42L, n_neighbors = 15L, min_dist = 0.3) {
  if (is.null(scores) || nrow(scores) < 4L) return(NULL)
  set.seed(seed)
  uwot::umap(scores, n_neighbors = min(n_neighbors, nrow(scores) - 1L),
             min_dist = min_dist, metric = "euclidean", verbose = FALSE)
}

## Wide data.frame wrapper: splits meta columns, transposes to features x cells,
## runs run_embedding(), returns data.frame of meta + UMAP1 + UMAP2.
run_umap_wide <- function(wide_df, meta_cols,
                          n_hvf = 1000L, n_pcs = 50L,
                          n_neighbors = 15L, seed = 42L) {
  if (is.null(wide_df) || nrow(wide_df) < 4L) return(NULL)
  meta <- wide_df %>% dplyr::select(dplyr::all_of(meta_cols))
  mat <- wide_df %>% dplyr::select(-dplyr::all_of(meta_cols)) %>% as.matrix()
  res <- run_embedding(t(mat), n_hvf = n_hvf, n_pcs = n_pcs,
                       n_neighbors = n_neighbors, seed = seed)
  if (is.null(res)) return(NULL)
  dplyr::bind_cols(meta, setNames(as.data.frame(res$umap), c("UMAP1", "UMAP2")))
}

## Per-column regression of score_mat on cell mean methylation.
## Both inputs are cells x features matrices. Returns residual matrix (cells x shared features).
regress_out_meth <- function(score_mat, meth_mat) {
  stopifnot(nrow(score_mat) == nrow(meth_mat))
  shared <- intersect(colnames(score_mat), colnames(meth_mat))
  if (length(shared) == 0) stop("regress_out_meth: no shared columns between score and meth matrices")
  score_mat <- score_mat[, shared, drop = FALSE]
  meth_mat <- meth_mat[, shared, drop = FALSE]
  cell_mean_meth <- rowMeans(meth_mat, na.rm = TRUE)
  resid_mat <- vapply(seq_len(ncol(score_mat)), function(j) {
    y <- score_mat[, j]
    fit <- try(lm(y ~ cell_mean_meth, na.action = na.exclude), silent = TRUE)
    if (inherits(fit, "try-error")) return(y)
    as.numeric(residuals(fit))
  }, numeric(nrow(score_mat)))
  colnames(resid_mat) <- shared
  resid_mat
}

## Mean silhouette width for cluster labels on a coordinate matrix.
## Returns NA when fewer than 2 unique labels are present.
sil_score <- function(coords, labels) {
  lv <- as.integer(as.factor(labels))
  if (length(unique(lv)) < 2L) return(NA_real_)
  mean(cluster::silhouette(lv, dist(coords))[, "sil_width"])
}
