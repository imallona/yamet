## Local test render for crc_embeddings.Rmd.
## Requires the corrected SCE rsynced from barbara:
##   data/crc/results/sce_windows_colon_corrected_10000.rds
##
## Run from the workflow/ directory:
##   cd ~/src/yamet/workflow && Rscript render_crc_embeddings_local.R

sce_path <- "data/crc/results/sce_windows_colon_corrected_10000.rds"
if (!file.exists(sce_path)) stop("corrected SCE not found at: ", sce_path)

snakemake <- readRDS("snmk_crc_embeddings.rds")
snakemake@input$feature_outputs <- list()
snakemake@threads <- 2L
snakemake@log$log <- "logs/render_crc_embeddings_local.log"

dir.create("logs", showWarnings = FALSE)

rmarkdown::render(
  "rules/src/crc_embeddings.Rmd",
  output_file   = normalizePath("crc_embeddings_local_test.html", mustWork = FALSE),
  knit_root_dir = normalizePath(".")
)
