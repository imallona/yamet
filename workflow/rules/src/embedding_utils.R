## Shared embedding utilities used by ecker.Rmd, argelaguet.Rmd, and
## crc_embeddings.Rmd. Source after ggtheme.R.
##
## Convention for run_embedding / run_pca_mat: mat is features x cells.
## Convention for run_umap_wide / regress_out_meth: cells x features (wide df).
## Requires dplyr, uwot, matrixStats, cluster.

## HVF selection, feature scaling, PCA (via randomized SVD), UMAP.
##
## NA handling strategy (two stages):
##
## Stage 1 (preferred): restrict HVF candidates to features with complete
## coverage across all cells (0 NAs). This keeps all cells and avoids
## imputation. NAs arise from missing sequencing coverage, not from the
## biological zero (unmethylated / zero entropy), which is a valid value.
## HVF score = var * coverage_frac; for 0-NA features coverage_frac = 1 so
## this reduces to plain variance.
## Stage 1 is used when at least min_complete_hvf complete features exist.
##
## Stage 2 (sparse fallback for feature sets like gene bodies):
## when too few complete features are available, HVFs are selected from the
## broader pool (features with <= max_na_frac missing cells) using
## coverage-weighted variance (score = var * coverage_frac), which penalises
## features that appear highly variable only because of sparse coverage.
## Cells missing at any selected HVF are then dropped rather than imputed.
## kept_cols records which cells were retained.
##
## Returns list(umap, pca, hvf_idx, kept_cols) or NULL.
## kept_cols is always returned: all TRUE in stage 1, subset in stage 2
## or after cell-level QC.
run_embedding <- function(mat,
                          n_hvf = 1000L,
                          n_pcs = 50L,
                          n_neighbors = 15L,
                          min_dist = 0.3,
                          seed = 42L,
                          max_na_frac = 0.3,
                          min_complete_hvf = 10L,
                          max_cell_na_frac = 0.5) {
  if (is.null(mat) || ncol(mat) < 4L) return(NULL)

  storage.mode(mat) <- "double"
  mat[!is.finite(mat)] <- NA

  ## Cell-level QC: drop cells with too many NAs before HVF selection.
  ## Analogous to scRNA-seq cell QC; cells with most loci unmeasured carry no
  ## useful signal and would otherwise inflate variance estimates.
  kept_qc <- colMeans(is.na(mat)) <= max_cell_na_frac
  if (sum(kept_qc) < 4L) {
    message("run_embedding: too few cells with <= ", max_cell_na_frac,
            " NA fraction; returning NULL")
    return(NULL)
  }
  if (sum(!kept_qc) > 0L)
    message("run_embedding: cell QC removed ", sum(!kept_qc),
            " cells with > ", max_cell_na_frac, " NA fraction")
  mat <- mat[, kept_qc, drop = FALSE]

  complete_rows <- rowSums(is.na(mat)) == 0L
  complete_mat <- mat[complete_rows, , drop = FALSE]
  vars_c <- matrixStats::rowVars(complete_mat, useNames = FALSE, na.rm = FALSE)
  n_complete <- sum(!is.na(vars_c) & vars_c > 0)

  if (n_complete >= as.integer(min_complete_hvf)) {
    rank_c <- order(vars_c, decreasing = TRUE)[seq_len(min(as.integer(n_hvf), n_complete))]
    hvf_idx <- which(complete_rows)[rank_c]
    sub_mat <- mat[hvf_idx, , drop = FALSE]
    kept_in_qc <- rep(TRUE, ncol(mat))
    message("run_embedding: stage 1 (complete features), HVF = ", length(hvf_idx),
            ", cells = ", ncol(mat))
  } else {
    message("run_embedding: only ", n_complete, " complete features",
            " (< min_complete_hvf = ", min_complete_hvf, ");",
            " sparse fallback (max_na_frac = ", max_na_frac, ")")
    pool_rows <- rowMeans(is.na(mat)) <= max_na_frac
    pool <- mat[pool_rows, , drop = FALSE]
    vars_p <- matrixStats::rowVars(pool, useNames = FALSE, na.rm = TRUE)
    cov_p <- rowMeans(!is.na(pool))
    score_p <- vars_p * cov_p
    n_keep <- min(as.integer(n_hvf), sum(!is.na(score_p) & score_p > 0))
    if (n_keep < 2L) {
      message("run_embedding: too few variable features in sparse fallback; returning NULL")
      return(NULL)
    }
    rank_p <- order(score_p, decreasing = TRUE)[seq_len(n_keep)]
    hvf_idx <- which(pool_rows)[rank_p]
    sub_mat <- mat[hvf_idx, , drop = FALSE]
    kept_in_qc <- colSums(is.na(sub_mat)) == 0L
    sub_mat <- sub_mat[, kept_in_qc, drop = FALSE]
    if (ncol(sub_mat) < 4L) {
      message("run_embedding: too few cells after NA removal in sparse fallback; returning NULL")
      return(NULL)
    }
    message("run_embedding: sparse fallback, HVF = ", length(hvf_idx),
            ", cells retained = ", sum(kept_in_qc), " / ", ncol(mat))
  }

  ## Map retained cells back to original column indices.
  kept_cols <- kept_qc
  kept_cols[kept_qc] <- kept_in_qc

  scaled <- t(scale(t(sub_mat)))
  bad <- apply(scaled, 1, function(x) all(is.na(x) | is.nan(x)))
  scaled <- scaled[!bad, , drop = FALSE]
  if (nrow(scaled) == 0L) return(NULL)

  n_pcs_use <- min(as.integer(n_pcs), nrow(scaled) - 1L, ncol(scaled) - 1L)
  if (n_pcs_use < 2L) {
    message("run_embedding: too few PCs available (", n_pcs_use, "); returning NULL")
    return(NULL)
  }
  set.seed(seed)
  pca_scores <- irlba::prcomp_irlba(t(scaled), n = n_pcs_use,
                                     center = FALSE, scale. = FALSE)$x

  nn <- min(as.integer(n_neighbors), nrow(pca_scores) - 1L)
  set.seed(seed)
  coords <- uwot::umap(pca_scores, n_neighbors = nn, min_dist = min_dist,
                       n_epochs = 500, verbose = FALSE)

  list(umap = coords, pca = pca_scores, hvf_idx = hvf_idx, kept_cols = kept_cols)
}

## PCA on a features x cells matrix. Returns cells x PCs score matrix or NULL.
run_pca_mat <- function(mat, n_pcs = 10L) {
  mat[!is.finite(mat)] <- NA
  mat <- mat[rowSums(is.na(mat)) == 0L, , drop = FALSE]
  if (nrow(mat) < 2L || ncol(mat) < 2L) return(NULL)
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
  if (is.null(scores) || nrow(scores) < 4L || ncol(scores) < 1L) return(NULL)
  if (all(!is.finite(scores))) return(NULL)
  set.seed(seed)
  uwot::umap(scores, n_neighbors = min(n_neighbors, nrow(scores) - 1L),
             min_dist = min_dist, metric = "euclidean", verbose = FALSE)
}

## Wide data.frame wrapper: splits meta columns, transposes to features x cells,
## runs run_embedding(), returns data.frame of meta + UMAP1 + UMAP2.
## In the sparse fallback, rows lacking HVF coverage are dropped; callers
## should note that the output may have fewer rows than the input.
run_umap_wide <- function(wide_df, meta_cols,
                          n_hvf = 1000L, n_pcs = 50L,
                          n_neighbors = 15L, seed = 42L,
                          max_na_frac = 0.3, min_complete_hvf = 10L,
                          max_cell_na_frac = 0.5) {
  if (is.null(wide_df) || nrow(wide_df) < 4L) return(NULL)
  meta <- wide_df %>% dplyr::select(dplyr::all_of(meta_cols))
  mat <- wide_df %>% dplyr::select(-dplyr::all_of(meta_cols)) %>% as.matrix()
  res <- run_embedding(t(mat), n_hvf = n_hvf, n_pcs = n_pcs,
                       n_neighbors = n_neighbors, seed = seed,
                       max_na_frac = max_na_frac,
                       min_complete_hvf = min_complete_hvf,
                       max_cell_na_frac = max_cell_na_frac)
  if (is.null(res)) return(NULL)
  dplyr::bind_cols(meta[res$kept_cols, , drop = FALSE],
                   setNames(as.data.frame(res$umap), c("UMAP1", "UMAP2")))
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
