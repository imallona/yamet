#!/usr/bin/env Rscript
##
## Parse SCNAs from CRC 2018 paper
##
## Atreya Choudhury
## Started 25th Sep 2025

library(readxl)

sheets <- excel_sheets(snakemake@input[[1]])
df_list <- lapply(
  sheets, \(sh) {
    excel <- read_excel(snakemake@input[[1]], sheet = sh, skip = 2, n_max = 3)
    df <- as.data.frame(t(excel[-1]), stringsAsFactors = FALSE, check.names = F)
    colnames(df) <- excel[[1]]
    df$cytoband <- rownames(df)
    rownames(df) <- NULL
    df$status <- strsplit(sh, "\\s+")[[1]][2]
    df
  }
)

scnas <- do.call(rbind, df_list)

bed_cols <- do.call(rbind, strsplit(scnas$`wide peak boundaries`, "[:-]"))
scnas$`wide peak boundaries` <- NULL

# assign as columns and convert start/end to numeric
scnas$chr <- bed_cols[, 1]
scnas$start <- as.numeric(bed_cols[, 2])
scnas$end <- as.numeric(bed_cols[, 3])

# move BED3 columns to the front
scnas <- scnas[
  , c("chr", "start", "end", setdiff(names(scnas), c("chr", "start", "end")))
]

# filter for SCNAs with q value < 0.05
scnas <- scnas[scnas$`residual q value` < 0.05, ]

tmp_bed <- tempfile(fileext = ".bed")
write.table(
  scnas, tmp_bed,
  col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE
)

system(
  paste0(
    "sort -k1,1 -k2,2n ", tmp_bed, " | gzip > ", snakemake@output[[1]]
  )
)

file.remove(tmp_bed)
