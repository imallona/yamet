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
#' @date 2025-04-21

options(scipen = 999)

samp <- as.integer(snakemake@wildcards[["sample"]])
f <- as.integer(snakemake@wildcards[["f"]])

set.seed(42 + samp)

stopifnot(f %% 8 == 1)

data_dir <- snakemake@params[["data"]]

template <- read.table(snakemake@input[[1]],
  sep = "\t",
  col.names = c("chr", "start", "end", "het", "snip_pos"),
  stringsAsFactors = F
)
template$snip_pos <- lapply(strsplit(template$snip_pos, ";"), as.numeric)

chain_gen <- function(n, snip_pos) {
  snips.length <- (n - 1) / 8
  chain <- c(rep("0011", snips.length), rep("0110", snips.length))
  # shuffle snips according to snip_pos
  chain[snip_pos] <- sample(chain[snip_pos])
  # flatten out snips sequence
  chain <- as.integer(unlist(strsplit(chain, split = "")))
  chain <- c(chain, 0)

  # set 30% of positions in a region to 0
  flip_num <- ceiling(30 * length(chain) / 100)
  flip_indices <- sample(seq_len(length(chain)), flip_num)
  chain[flip_indices] <- 0

  return(chain)
}

result <- do.call(rbind, lapply(seq_len(nrow(template)), function(i) {
  row <- as.list(template[i, ])
  chain <- data.frame(
    chr = rep(row$chr, row$end - row$start),
    pos = seq(row$start, by = 1, length.out = row$end - row$start),
    total = rep(1, row$end - row$start)
  )
  chain$beta <- chain_gen(row$end - row$start, row$snip_pos[[1]])
  chain$meth <- chain$beta
  return(chain)
}))
write.table(
  result[, c("chr", "pos", "meth", "total", "beta")],
  snakemake@output[[1]],
  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)
