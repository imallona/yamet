suppressPackageStartupMessages({
  library(markovchain)
})

options(scipen = 999)
set.seed(42)

samples <- as.integer(snakemake@wildcards[["samples"]])
data_dir <- snakemake@params[["data"]]

template <- read.table(snakemake@input[[1]],
  sep = "\t",
  col.names = c("chr", "start", "end", "rho", "prob.0"),
)
N <- nrow(template)
f <- template$end[1] - template$start[1]

mcgen <- function(n, rho, prob.0) {
  tm <- matrix(c(
    rho, 1 - rho,
    1 - rho, rho
  ), nrow = 2, byrow = T)

  mc <- new("markovchain",
    states = c("0", "1"),
    transitionMatrix = tm
  )
  return(as.integer(
    rmarkovchain(n, mc, t0 = sample(c("0", "1"), 1, replace = T, prob = c(prob.0, 1 - prob.0)))
  ))
}

genSamp <- function(sample, template) {
  result <- do.call(rbind, lapply(seq_len(nrow(template)), function(i) {
    row <- as.list(template[i, ])
    chain <- data.frame(
      chr = rep(row$chr, row$end - row$start),
      pos = seq(row$start, by = 1, length.out = row$end - row$start),
      total = rep(1, row$end - row$start)
    )
    chain$beta <- mcgen(row$end - row$start, row$rho, row$prob.0)
    chain$meth <- chain$beta
    return(chain)
  }))
  write.table(
    result[, c("chr", "pos", "meth", "total", "beta")],
    paste0(data_dir, paste("/sim", sample, samples, N, f, "tsv", sep = ".")),
    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
  )
}

for (sample in seq_len(samples)) {
  genSamp(sample, template)
}
