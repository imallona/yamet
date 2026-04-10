## Shared helpers for building SingleCellExperiment objects from yamet
## windows outputs (det.out.gz, norm.det.out.gz, meth.out.gz).
##
## Used by argelaguet_windows.Rmd and ecker_windows.Rmd.
## Requires: readr, SingleCellExperiment, SummarizedExperiment, BiocParallel.
##
## Column layout in yamet windows output files:
##   chr, start, end, <input_cell_path_1>, <input_cell_path_2>, ...
## -1 encodes missing coverage (read as NA via na = "-1").

## Read one windows output file (det, norm.det, or meth).
## Returns list(coords, mat) where mat is regions x cells.
## cell_id_fun transforms full column names (file paths) to cell IDs.
read_yamet_det <- function(filepath, na_val = "-1",
                           cell_id_fun = function(x) sub("\\.gz$", "", basename(x))) {
  dt <- readr::read_tsv(filepath, na = na_val, show_col_types = FALSE)
  coord_cols <- c("chr", "start", "end")
  cell_cols <- setdiff(colnames(dt), coord_cols)
  coords <- dt[, coord_cols]
  coords$region <- sprintf("%s:%s-%s", coords$chr, coords$start, coords$end)
  mat <- as.matrix(dt[, cell_cols, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- coords$region
  colnames(mat) <- cell_id_fun(cell_cols)
  list(coords = coords, mat = mat)
}

## Build a SingleCellExperiment from a named list of windows output file paths.
##
## assay_paths: named list mapping assay names to file paths, e.g.
##   list(sampen = "...det.out.gz", meth = "...meth.out.gz",
##        adjS   = "...norm.det.out.gz")
##
## col_data: data.frame with rownames = cell IDs. Cells not present in ALL
##   assay files are dropped without imputation.
build_windows_sce <- function(assay_paths, col_data,
                              cell_id_fun = function(x) sub("\\.gz$", "", basename(x))) {
  dets <- lapply(assay_paths, read_yamet_det, cell_id_fun = cell_id_fun)
  shared_cells <- Reduce(intersect,
                         c(list(rownames(col_data)),
                           lapply(dets, function(x) colnames(x$mat))))
  shared_regions <- Reduce(intersect, lapply(dets, function(x) rownames(x$mat)))
  if (length(shared_cells) < 2L)
    stop("build_windows_sce: fewer than 2 cells shared across all assay files and col_data")
  assay_list <- lapply(dets, function(x)
    x$mat[shared_regions, shared_cells, drop = FALSE])
  SingleCellExperiment::SingleCellExperiment(
    assays  = assay_list,
    colData = col_data[shared_cells, , drop = FALSE]
  )
}

## Fit per-window sampen ~ meth + I(meth^2) and store residuals as
## sampen_corrected. NAs in sampen or meth are passed through to the residuals.
## BPPARAM: BiocParallel parameter for parallelism.
add_meth_correction <- function(sce, BPPARAM = BiocParallel::SerialParam()) {
  rowwise_residuals <- function(i) {
    df <- data.frame(sampen = SummarizedExperiment::assay(sce, "sampen")[i, ],
                     meth   = SummarizedExperiment::assay(sce, "meth")[i, ])
    fit <- try(lm(sampen ~ meth + I(meth^2), data = df), silent = TRUE)
    if (inherits(fit, "try-error")) return(rep(NA_real_, ncol(sce)))
    as.numeric(residuals(fit))
  }
  res_mat <- do.call(rbind, BiocParallel::bplapply(
    seq_len(nrow(sce)), rowwise_residuals, BPPARAM = BPPARAM))
  rownames(res_mat) <- rownames(sce)
  colnames(res_mat) <- colnames(sce)
  SummarizedExperiment::assay(sce, "sampen_corrected") <- res_mat
  sce
}
