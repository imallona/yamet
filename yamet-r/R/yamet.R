#' Main yamet analysis function
#'
#' Reads methylation data from multiple input files, computes sample entropies,
#' and returns a SummarizedExperiment object containing the results.
#'
#' @param filenames Character vector of input file paths.
#'   Each file must be tab-separated, sorted by chromosome and position, and
#'   contain five columns: chromosome, position, methylated reads, total reads and rate.
#'   Optionally, files can be gzipped.
#' @param reference_path Path to the reference file.
#'   The file should be tab-separated, sorted by chromosome and position, and
#'   contain two columns: chromosome and position.
#'   Optionally, the file can be gzipped.
#' @param intervals_path Path to BED file or similar defining genomic intervals.
#'   It should be tab-separated, sorted by chromosome and start position, and
#'   contain three columns: chromosome, start and end.
#' @param cores Number of cores to use (default: 1).
#' @param chunk_size per file chunk size bytes as an integer (default: 64K).
#'   Number of bytes to process at once in each file.
#' @param skip_header_cell number of lines skipped in cell files (default: 0).
#' @param skip_header_reference number of lines skipped in the reference file (default: 0).
#' @param skip_header_intervals number of lines skipped in the intervals file (default: 0).
#'
#' @return A \link[SummarizedExperiment]{SummarizedExperiment} object
#'   containing assays (`sampens`, `meths`), row metadata (genomic intervals),
#'   and sample-level colData.
#'
#' @export
yamet <- function(
    filenames = c(), reference_path = NULL, intervals_path = NULL, cores = NULL,
    chunk_size = 64 * 1024, skip_header_cell = 0,
    skip_header_reference = 0, skip_header_intervals = 0) {
  if (length(filenames) == 0 || !all(file.exists(filenames))) {
    stop("Please provide at least one cell file in 'filenames'.")
  }
  if (length(reference_path) != 1 || !file.exists(reference_path)) {
    stop("Please provide a reference file in 'reference_path'.")
  }
  if (length(intervals_path) != 1 || !file.exists(intervals_path)) {
    stop("Please provide a intervals file in 'intervals_path'.")
  }

  max_cores <- parallel::detectCores(logical = FALSE)
  if (missing(cores) || is.null(cores)) {
    cores <- max(1, floor(max_cores - log2(max_cores)))
  } else {
    cores <- min(as.integer(cores), max_cores)
  }

  t <- system.time({
    yamet_raw <- yamet_cpp(
      filenames,
      reference_path,
      intervals_path,
      cores,
      as.integer(chunk_size),
      skip_header_cell,
      skip_header_reference,
      skip_header_intervals
    )
  })

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      sampens = yamet_raw$sampens,
      meths = yamet_raw$meths
    ),
    rowData = data.frame(
      chr = yamet_raw$chr,
      start = yamet_raw$start,
      end = yamet_raw$end,
      shannon = yamet_raw$shannon,
      meth = yamet_raw$avg_meth
    ),
    colData = data.frame(
      sample = filenames,
      sampen = yamet_raw$cell_sampen,
      meth = yamet_raw$cell_meth
    ),
    metadata = list(
      reference = reference_path,
      intervals = intervals_path,
      processing_date = Sys.Date(),
      t_user = t["user.self"],
      t_system = t["sys.self"],
      t_elapsed = t["elapsed"]
    )
  )
  rownames(se) <- paste0(
    yamet_raw$chr, ":", yamet_raw$start, "-", yamet_raw$end
  )
  colnames(se) <- filenames
  return(se)
}
