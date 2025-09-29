#!/usr/bin/env Rscript
##
## Parse SCNAs from patient 1 of
# Single-cell multiomics sequencing and analyses of human colorectal cancer

# Shuhui Bian https://orcid.org/0000-0002-9662-113X, Yu Hou https://orcid.org/0000-0001-7875-7087, Xin Zhou https://orcid.org/0000-0002-4048-4017, Xianlong Li https://orcid.org/0000-0001-7745-2695, Jun Yong https://orcid.org/0000-0002-3770-2108, Yicheng Wang https://orcid.org/0000-0002-3901-1482, Wendong Wang https://orcid.org/0000-0002-8631-3109, Jia Yan https://orcid.org/0000-0002-6203-0016, Boqiang Hu https://orcid.org/0000-0002-2045-3619, Hongshan Guo https://orcid.org/0000-0003-2799-2989, Jilian Wang https://orcid.org/0000-0002-3951-6931, Shuai Gao https://orcid.org/0000-0002-8208-9197, Yunuo Mao https://orcid.org/0000-0003-1428-270X, Ji Dong, Ping Zhu, Dianrong Xiu, Liying Yan https://orcid.org/0000-0001-9572-9440, Lu Wen, Jie Qiao https://orcid.org/0000-0003-2126-1376 , Fuchou Tang https://orcid.org/0000-0002-8625-7717 , and Wei Fu https://orcid.org/0000-0001-5248-7891


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
