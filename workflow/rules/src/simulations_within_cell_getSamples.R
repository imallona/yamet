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

s <- as.integer(snakemake@wildcards[["s"]])
f <- as.integer(snakemake@wildcards[["f"]])

set.seed(42 + s)

stopifnot(f %% 8 == 1)

data_dir <- snakemake@params[["data"]]

template <- read.table(snakemake@input[[1]],
  header = T, stringsAsFactors = F
)
template$snip_pos <- lapply(strsplit(template$snip_pos, ";"), as.integer)
template$delta <- lapply(strsplit(template$delta, ";"), as.integer)

chain_gen <- function(n, snip_pos, higher, delta) {
  snips.length <- (n - 1) / 8
  chain <- c(rep("0011", snips.length), rep("0110", snips.length))
  # shuffle snips according to snip_pos
  chain[snip_pos] <- sample(chain[snip_pos])
  # sample 35% of snips
  snips_subset_count <- ceiling(0.35 * length(chain))
  snips_subset_idx <- sample(seq_len(length(chain)), snips_subset_count)
  # select snips to methylate/demethylate with 'delta' differential
  lower_flip_count <- floor(snips_subset_count * (1 - delta * 0.09) / 2)
  flip_indices <- sample(snips_subset_idx, lower_flip_count)
  # methylate/demethylate according to 'higher'
  chain[flip_indices] <- paste0(c(rep(1 - higher, 3), higher), collapse = "")
  chain[setdiff(snips_subset_idx, flip_indices)] <- paste0(
    c(rep(higher, 3), 1 - higher),
    collapse = ""
  )
  # flatten out snips sequence
  chain <- as.integer(unlist(strsplit(chain, split = "")))
  chain <- c(chain, 0)

  return(chain)
}

result <- do.call(rbind, lapply(seq_len(nrow(template)), function(i) {
  row <- as.list(template[i, ])
  chain <- data.frame(
    chr = rep(row$chr, row$end - row$start),
    pos = seq(row$start, by = 1, length.out = row$end - row$start),
    total = rep(1, row$end - row$start)
  )
  chain$beta <- chain_gen(
    row$end - row$start, row$snip_pos[[1]],
    row$higher, row$delta[[1]][s]
  )
  chain$meth <- chain$beta
  return(chain)
}))
write.table(
  result[, c("chr", "pos", "meth", "total", "beta")],
  snakemake@output[[1]],
  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)
