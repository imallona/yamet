#' Simulates a template for genomic feature regions
#'
#' Inputs (via Snakemake wildcards):
#'   - S: Number of samples
#'   - N: Number of feature regions
#'   - f: Length of each feature region (must satisfy f %% 8 == 1)
#'
#' It generates two tab-delimited output files:
#'   - A reference file listing all positions across all features.
#'   - An annotation file with chromosome, region bounds,
#'     a sampled variability index, and indices used for shuffling
#'     snips within features (described in the strategy)
#'
#' @author Atreya Choudhury
#' @date 2025-04-21

options(scipen = 999)
set.seed(42)

S <- as.integer(snakemake@wildcards[["S"]])
N <- as.integer(snakemake@wildcards[["N"]])
f <- as.integer(snakemake@wildcards[["f"]])

stopifnot(f %% 8 == 1)

out <- snakemake@output[[1]]
out_dir <- dirname(snakemake@output[[1]])
dir.create(out_dir, recursive = TRUE)

# reference generation

total.positions <- N * f

chr <- rep("chrSim", total.positions)
pos <- seq(0, total.positions - 1)

write.table(
  data.frame(chr = chr, pos = pos),
  paste(out_dir, "/ref.", N, ".", f, ".tsv", sep = ""),
  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)

# template generation

chr <- rep("chrSim", N)
start <- seq(0, by = f, length.out = N)
end <- start + f
# Assign each region a vi score
vi <- sample(1:10, size = N, replace = TRUE)

# The number of 0011 or 0110 snips in every feature
snips.length <- (f - 1) / 8

# positions of snips we will shuffle
snip_pos <- sapply(vi, function(x) {
  # fraction of 0011 and 0110 snips to shuffle
  shuf_count <- floor(x * snips.length / (1.5 * 10))
  # select an equal fraction given by shuf_count from 0011 and 0110 snips
  # and record their snip indices
  paste(sort(c(
    sample(1:snips.length, shuf_count, replace = F),
    snips.length + sample(1:snips.length, shuf_count, replace = F)
  )), collapse = ";")
})

# variable deciding whether in a sample, more 0s are flipped or 1s
higher <- sample(0:1, size = N, replace = TRUE)

# variable quantifying by how much more 0s/1s are flipped in a sample
delta <- sapply(vi, function(x) {
  paste(
    extraDistr::rbbinom(S, 10, alpha = 11 - x, beta = 8),
    collapse = ";"
  )
})

write.table(
  data.frame(
    chr = chr, start = start, end = end,
    vi = vi, snip_pos = snip_pos,
    higher = higher, delta = delta
  ),
  out,
  sep = "\t", row.names = FALSE, quote = FALSE
)
