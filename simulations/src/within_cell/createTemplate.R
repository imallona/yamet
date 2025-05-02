options(scipen = 999)
set.seed(42)

features.N <- as.integer(snakemake@wildcards[["N"]])
features.length <- as.integer(snakemake@wildcards[["f"]])

stopifnot(features.length %% 8 == 1)

out <- snakemake@output[[1]]
out_dir <- dirname(snakemake@output[[1]])
dir.create(out_dir, recursive = TRUE)

total.positions <- features.N * features.length

chr <- rep("chrSim", total.positions)
pos <- seq(0, total.positions - 1)

write.table(
  data.frame(chr = chr, pos = pos),
  paste(out_dir, "/ref.", features.N, ".", features.length, ".tsv", sep = ""),
  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)

chr <- rep("chrSim", features.N)
lbound <- seq(0, by = features.length, length.out = features.N)
ubound <- lbound + features.length
vi <- sample(1:10, size = features.N, replace = TRUE)

snips.length <- (features.length - 1) / 8

snip_pos <- sapply(vi, function(x) {
  shuf_count <- floor(x * snips.length / (1.5 * 10))
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
