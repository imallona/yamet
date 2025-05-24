#' Main yamet analysis function
#'
#' Reads methylation data from multiple input files, computes sample entropies,
#' and returns a SummarizedExperiment object containing the results.
#'
#' @param filenames Character vector of input files.
#' @param reference_path Path to reference data file.
#' @param intervals_path Path to BED file or similar defining genomic intervals.
#' @param cores Number of cores to use (default: 1).
#' @param chunk_size Processing chunk size (default: 1000).
#'
#' @return A \link[SummarizedExperiment]{SummarizedExperiment} object
#'   containing assays (`sampens`, `meths`), row metadata (genomic intervals),
#'   and sample-level colData.
#'
#' @export


yamet <- function(
    filenames, reference_path, intervals_path, cores = 1, chunk_size = 1000) {
  if (!all(file.exists(filenames))) {
    stop("Some input files don't exist")
  }

  yamet_raw <- yamet_cpp(
    filenames,
    reference_path,
    intervals_path,
    as.integer(cores),
    as.integer(chunk_size)
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      sampens = yamet_raw$sampens,
      meths = yamet_raw$meths
    ),
    rowData = data.frame(
      chr = yamet_raw$chr,
      start = yamet_raw$start,
      end = yamet_raw$end
    ),
    colData = data.frame(
      sample = filenames,
      sampen = yamet_raw$avg_sampen,
      meth = yamet_raw$avg_meth
    ),
    metadata = list(
      reference = reference_path,
      intervals = intervals_path,
      processing_date = Sys.Date()
    )
  )
  rownames(se) <- paste0(
    yamet_raw$chr, ":", yamet_raw$start, "-", yamet_raw$end
  )
  colnames(se) <- filenames
  return(se)
}
