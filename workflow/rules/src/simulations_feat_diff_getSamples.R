suppressPackageStartupMessages({
  library(markovchain)
})

options(scipen = 999)

samples <- as.integer(snakemake@params[["samples"]])
data_dir <- snakemake@params[["data"]]

template <- read.table(snakemake@input[[1]],
  sep = "\t",
  col.names = c("chr", "start", "end", "type", "prob.0"),
)
N <- nrow(template)
f <- template$end[1] - template$start[1]

set.seed(42)

mcgen <- function(n, type, var, prob.0) {
  if (type == "lmr") {
    base <- 0.8
    tm <- matrix(c(
      base, 1 - base,
      base, 1 - base
    ), nrow = 2, byrow = T)
  } else if (type == "hmr") {
    base <- 0.2
    tm <- matrix(c(
      base, 1 - base,
      base, 1 - base
    ), nrow = 2, byrow = T)
  } else if (type == "imrCons") {
    base <- 0.8
    tm <- matrix(c(
      base, 1 - base,
      1 - base, base
    ), nrow = 2, byrow = T)
  } else if (type == "imrRand") {
    base <- 0.5
    tm <- matrix(c(
      base, 1 - base,
      1 - base, base
    ), nrow = 2, byrow = T)
  }

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
    chain$beta <- mcgen(row$end - row$start, row$type, row$var, row$prob.0)
    chain$meth <- chain$beta
    return(chain)
  }))
  write.table(
    result[, c("chr", "pos", "meth", "total", "beta")],
    paste(data_dir, "/sim.", N, ".", f, ".", sample, ".tsv", sep = ""),
    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
  )
}

for (sample in seq_len(samples)) {
  genSamp(sample, template)
}
