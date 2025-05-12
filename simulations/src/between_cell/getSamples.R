#' Generates one sample using a template generated in createTemplate.R
#'
#' Inputs (via Snakemake wildcards):
#'   - N: Number of feature regions
#'   - f: Length of each feature region (must satisfy f %% 8 == 1)
#'
#' It generates a tab-delimited output file:
#'   - A cell sample file with chromosome, genomic position,
#'     methylated reads, total reads and methylation rate
#'
#' @author Atreya Choudhury
#' @date 2025-04-17

suppressPackageStartupMessages({
  library(markovchain)
})

options(scipen = 999)

s <- as.integer(snakemake@wildcards[["s"]])

set.seed(42 + s)

data_dir <- snakemake@params[["data"]]

template <- read.table(snakemake@input[[1]],
  header = T, stringsAsFactors = F
)
template$rho <- lapply(strsplit(template$rho, ";"), as.numeric)

mcgen <- function(n, rho) {
  tm <- matrix(c(
    rho, 1 - rho,
    rho, 1 - rho
  ), nrow = 2, byrow = T)

  mc <- new("markovchain",
    states = c("0", "1"),
    transitionMatrix = tm
  )
  return(as.integer(
    rmarkovchain(n, mc, t0 = sample(c("0", "1"), 1, replace = T))
  ))
}

result <- do.call(rbind, lapply(seq_len(nrow(template)), function(i) {
  row <- as.list(template[i, ])
  chain <- data.frame(
    chr = rep(row$chr, row$end - row$start),
    pos = seq(row$start, by = 1, length.out = row$end - row$start),
    total = rep(1, row$end - row$start)
  )
  chain$beta <- mcgen(row$end - row$start, row$rho[[1]][s])
  chain$meth <- chain$beta
  return(chain)
}))
write.table(
  result[, c("chr", "pos", "meth", "total", "beta")],
  snakemake@output[[1]],
  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)
