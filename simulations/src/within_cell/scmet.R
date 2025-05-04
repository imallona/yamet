suppressPackageStartupMessages({
  library(data.table)
})

samples <- as.integer(snakemake@wildcards[["samples"]])
N <- as.integer(snakemake@wildcards[["N"]])
f <- as.integer(snakemake@wildcards[["f"]])
data_dir <- dirname(snakemake@output[[1]])

intervals <- read.table(
  paste0(data_dir, paste("/intervals", N, f, "tsv", sep = ".")), ,
  col.names = c("chr2", "start2", "end2", "vi", "snip_pos"),
  header = F
)
meth <- read.table(
  paste0(data_dir, paste("/yamet", samples, N, f, "meth.out", sep = ".")),
  header = TRUE
)
jnt <- cbind(meth, intervals)
output_cols <- grep("^output", colnames(jnt), value = TRUE)
scmet_Y <- do.call(rbind, lapply(output_cols, function(col) {
  data.table(
    Feature = paste0(jnt$start, ":", jnt$end),
    Cell = col,
    total_reads = jnt$end - jnt$start,
    met_reads = round(jnt[[col]] * (jnt$end - jnt$start))
  )
}))
scmet_fit <- suppressWarnings(
  scMET::scmet(Y = scmet_Y, L = 4, iter = 2000, n_cores = snakemake@threads, seed = 42)
)
saveRDS(scmet_fit, snakemake@output[[1]])
