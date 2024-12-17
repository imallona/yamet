#!/usr/bin/env Rscript
##
## Harmonizes Ecker's 2021 nemo and paper metadata to get sane cell annotations
##
## Izaskun Mallona
## Started 4th Dec 2024

library(readxl)

# print(snakemake)

## xlsx thing
paper_fn <- snakemake@input$paper

excel <- as.data.frame(read_excel(paper_fn, skip = 14))
colnames(excel)[1] <- 'cell_nemo'

## tsv, gz compressed
nemo_fn <-  snakemake@input$nemo

nemo <- read.table(nemo_fn, header = TRUE)

merged <- merge(x = excel, y = nemo, by.x = 'cell_nemo', by.y = 'cell', all.y = TRUE) 

stopifnot(nrow(merged) == 9941)

merged$basename <- basename(merged$AllcPath)

## write, gz-compressed
gz <- gzfile(snakemake@output$metadata, "w")
write.table(merged, gz, col.names = TRUE, row.names = FALSE, sep = "\t")
close(gz)
