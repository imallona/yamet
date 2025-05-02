#' Simulates a template for genomic feature regions
#'
#' Inputs (via Snakemake wildcards):
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

features.N <- as.integer(snakemake@wildcards[["N"]])
features.length <- as.integer(snakemake@wildcards[["f"]])

stopifnot(features.length %% 8 == 1)

out <- snakemake@output[[1]]
out_dir <- dirname(snakemake@output[[1]])
dir.create(out_dir, recursive = TRUE)

# reference generation

total.positions <- features.N * features.length

chr <- rep("chrSim", total.positions)
pos <- seq(0, total.positions - 1)

write.table(
  data.frame(chr = chr, pos = pos),
  paste(out_dir, "/ref.", features.N, ".", features.length, ".tsv", sep = ""),
  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)

# template generation

chr <- rep("chrSim", features.N)
lbound <- seq(0, by = features.length, length.out = features.N)
ubound <- lbound + features.length
# Assign each region a vi score
vi <- sample(1:10, size = features.N, replace = TRUE)

# The number of 0011 or 0110 snips in every feature
snips.length <- (features.length - 1) / 8

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

write.table(
  data.frame(
    chr = chr, lbound = lbound, ubound = ubound,
    vi = vi, snip_pos = snip_pos
  ),
  out,
  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)
