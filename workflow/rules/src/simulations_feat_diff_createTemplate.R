options(scipen = 999)

features.N <- as.integer(snakemake@wildcards[["N"]])
features.length <- as.integer(snakemake@wildcards[["f"]])
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

types <- c("lmr", "hmr", "imrCons", "imrRand")

chr <- rep("chrSim", features.N)
lbound <- seq(0, by = features.length, length.out = features.N)
ubound <- lbound + features.length
features <- sample(types, size = features.N, replace = TRUE)
prob.0 <- runif(features.N, 0, 1)

write.table(
  data.frame(
    chr = chr, lbound = lbound, ubound = ubound,
    features = features, prob.0 = prob.0
  ),
  out,
  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)
