## Shared embedding utilities used by ecker.Rmd, argelaguet.Rmd, and
## crc_embeddings.Rmd. Source after ggtheme.R.
##
## Convention for run_embedding / run_pca_mat: mat is features x cells.
## Convention for run_umap_wide / regress_out_meth: cells x features (wide df).
## Requires dplyr, uwot, matrixStats, cluster, pcaMethods.

## HVF selection, feature scaling, PCA, UMAP. NA-aware throughout.
##
## scBS-seq methylation and entropy values have a meaningful zero (unmethylated,
## or zero local entropy) which is distinct from missing coverage (NA). Standard
## scRNA-seq pipelines treat zero as missing, which is wrong here, so HVF /
## PCA / UMAP must keep NAs distinct from zeros and never impute.
##
## Two-stage strategy:
##
## Stage 1 (preferred): when at least min_complete_hvf features have complete
## coverage across all post-QC cells, HVF candidates are restricted to those
## features and ranked by plain variance. PCA is randomized SVD via irlba on
## the row-scaled matrix. All post-QC cells are kept.
##
## Stage 2 (sparse fallback): when too few complete features exist, HVFs are
## drawn from the broader pool of features with at most max_na_frac missing
## cells, ranked by coverage-weighted variance (score = var * coverage_frac)
## to penalise features whose apparent variance comes from sparse coverage.
## NAs in the selected HVF block are NOT imputed and cells are NOT dropped.
## PCA is NIPALS via pcaMethods, which handles missing values natively by
## iteratively regressing each component on the available entries (no
## explicit fill of missing values). UMAP runs on the resulting NIPALS
## scores, which are finite for every cell.
##
## Cells are only filtered by max_cell_na_frac (cell QC). Whichever stage is
## taken, kept_cols reports the original column positions surviving QC.
##
## Returns list(umap, pca, hvf_idx, kept_cols, stage) or NULL.
run_embedding <- function(mat,
                          n_hvf = 1000L,
                          n_pcs = 50L,
                          n_neighbors = 15L,
                          min_dist = 0.3,
                          seed = 42L,
                          max_na_frac = 0.3,
                          min_complete_hvf = 100L,
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

  use_nipals <- FALSE
  if (n_complete >= as.integer(min_complete_hvf)) {
    stage <- "stage1_complete"
    rank_c <- order(vars_c, decreasing = TRUE)[seq_len(min(as.integer(n_hvf), n_complete))]
    hvf_idx <- which(complete_rows)[rank_c]
    sub_mat <- mat[hvf_idx, , drop = FALSE]
    message("run_embedding: stage 1 (complete features), HVF = ", length(hvf_idx),
            ", cells = ", ncol(mat))
  } else {
    stage <- "stage2_nipals"
    use_nipals <- TRUE
    message("run_embedding: only ", n_complete, " complete features",
            " (< min_complete_hvf = ", min_complete_hvf, ");",
            " sparse fallback with NIPALS PCA (max_na_frac = ",
            max_na_frac, ")")
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
    message("run_embedding: sparse fallback, HVF = ", length(hvf_idx),
            ", cells = ", ncol(mat),
            " (NAs retained for NIPALS PCA, no cell dropping)")
  }

  ## Row-wise centering and scaling. na.rm preserves NAs in stage 2 so NIPALS
  ## sees the missingness pattern, while in stage 1 the rows have no NAs.
  row_mean <- matrixStats::rowMeans2(sub_mat, na.rm = TRUE)
  row_sd <- matrixStats::rowSds(sub_mat, na.rm = TRUE)
  bad_row <- !is.finite(row_sd) | row_sd == 0
  if (any(bad_row)) {
    sub_mat <- sub_mat[!bad_row, , drop = FALSE]
    hvf_idx <- hvf_idx[!bad_row]
    row_mean <- row_mean[!bad_row]
    row_sd <- row_sd[!bad_row]
  }
  if (nrow(sub_mat) < 2L) {
    message("run_embedding: no usable HVFs after row scaling; returning NULL")
    return(NULL)
  }
  scaled <- (sub_mat - row_mean) / row_sd

  n_pcs_use <- min(as.integer(n_pcs), nrow(scaled) - 1L, ncol(scaled) - 1L)
  if (n_pcs_use < 2L) {
    message("run_embedding: too few PCs available (", n_pcs_use, "); returning NULL")
    return(NULL)
  }

  set.seed(seed)
  if (use_nipals) {
    ## pcaMethods::pca with method = "nipals" handles missing values via
    ## iterative regression on observed entries; no imputation step replaces
    ## NAs with point estimates before factorisation. Input expected as
    ## samples (cells) x variables (features).
    pc_obj <- pcaMethods::pca(t(scaled), method = "nipals",
                              nPcs = n_pcs_use, scale = "none",
                              center = FALSE, seed = seed)
    pca_scores <- pcaMethods::scores(pc_obj)
    if (any(!is.finite(pca_scores))) {
      message("run_embedding: NIPALS produced non-finite scores; returning NULL")
      return(NULL)
    }
  } else if (nrow(scaled) <= 2L * n_pcs_use) {
    ## Few features relative to requested PCs: use base prcomp to avoid the
    ## irlba "too large a percentage of total singular values" warning.
    pca_scores <- prcomp(t(scaled), center = FALSE, scale. = FALSE,
                         rank. = n_pcs_use)$x
  } else {
    pca_scores <- irlba::prcomp_irlba(t(scaled), n = n_pcs_use,
                                       center = FALSE, scale. = FALSE)$x
  }

  nn <- min(as.integer(n_neighbors), nrow(pca_scores) - 1L)
  set.seed(seed)
  coords <- uwot::umap(pca_scores, n_neighbors = nn, min_dist = min_dist,
                       n_epochs = 500, verbose = FALSE)

  list(umap = coords, pca = pca_scores, hvf_idx = hvf_idx,
       kept_cols = kept_qc, stage = stage)
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
## Cells failing the max_cell_na_frac QC are dropped; in both the complete-feature
## branch and the NIPALS sparse fallback no further cells are removed.
run_umap_wide <- function(wide_df, meta_cols,
                          n_hvf = 1000L, n_pcs = 50L,
                          n_neighbors = 15L, seed = 42L,
                          max_na_frac = 0.3, min_complete_hvf = 100L,
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

## NA-aware per-feature variance explained by a grouping variable.
## mat: features x cells, groups: factor/character of length ncol(mat).
## Returns a named numeric vector (one R-squared per feature).
row_variance_explained <- function(mat, groups) {
  groups <- as.factor(groups)
  vapply(seq_len(nrow(mat)), function(i) {
    y <- mat[i, ]
    ok <- !is.na(y)
    if (sum(ok) < 3L || length(unique(groups[ok])) < 2L) return(NA_real_)
    ss_tot <- sum((y[ok] - mean(y[ok]))^2)
    if (ss_tot == 0) return(NA_real_)
    gm <- tapply(y[ok], groups[ok], mean)
    ss_between <- sum((gm[as.character(groups[ok])] - mean(y[ok]))^2)
    ss_between / ss_tot * 100
  }, numeric(1))
}

## Mean silhouette width for cluster labels on a coordinate matrix.
## Returns NA when fewer than 2 unique labels are present.
sil_score <- function(coords, labels) {
  lv <- as.integer(as.factor(labels))
  if (length(unique(lv)) < 2L) return(NA_real_)
  mean(cluster::silhouette(lv, dist(coords))[, "sil_width"])
}
