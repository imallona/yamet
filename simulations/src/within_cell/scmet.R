suppressPackageStartupMessages({
  library(data.table)
})

S <- as.integer(snakemake@wildcards[["S"]])
N <- as.integer(snakemake@wildcards[["N"]])
f <- as.integer(snakemake@wildcards[["f"]])
data_dir <- dirname(snakemake@output[[1]])

intervals <- fread(
  paste0(data_dir, paste("/intervals", S, N, f, "tsv", sep = ".")),
  header = T, select = "vi"
)
meth <- fread(
  paste0(data_dir, paste("/yamet", S, N, f, "meth.out", sep = ".")),
  header = T
)
jnt <- cbind(meth, intervals)
output_cols <- grep("^output", colnames(jnt), value = TRUE)
scmet_Y <- melt(jnt,
  measure.vars = output_cols,
  variable.name = "Cell",
  value.name = "meth_pct"
)[
  , .(
    Feature = paste0(start, ":", end),
    Cell,
    total_reads = end - start,
    met_reads = round(meth_pct * (end - start))
  )
]
scmet_fit <- suppressWarnings(
  scMET::scmet(Y = scmet_Y, L = 4, iter = 5000, n_cores = snakemake@threads, seed = 42)
)
saveRDS(scmet_fit, snakemake@output[[1]])
