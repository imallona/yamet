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
#' @date 2025-04-17

options(scipen = 999)
set.seed(42)

S <- as.integer(snakemake@wildcards[["S"]])
N <- as.integer(snakemake@wildcards[["N"]])
f <- as.integer(snakemake@wildcards[["f"]])

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

rho <- sapply(vi, function(x) {
  paste(
    sample(seq(0.1, 0.9, length.out = 10)[1:x], S, replace = TRUE),
    collapse = ";"
  )
})

write.table(
  data.frame(
    chr = chr, start = start, end = end,
    bc_vi = vi, rho = rho
  ),
  out,
  sep = "\t", row.names = FALSE, quote = FALSE
)
