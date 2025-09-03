suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(svglite)
})

bet_rds <- readRDS(snakemake@params[["bet_rds_path"]])
wit_rds <- readRDS(snakemake@params[["wit_rds_path"]])

arranged <- ggarrange(
  ggarrange(
    # top row: sampen vs meth and shannon
    ggarrange(
      bet_rds$main_plots$samp_meth,
      bet_rds$main_plots$samp_shan,
      wit_rds$main_plots$samp_meth,
      wit_rds$main_plots$samp_shan,
      legend = "none",
      nrow = 1
    ),
    ggarrange(
      ggarrange(
        bet_rds$main_plots$scmet_mu_meth,
        wit_rds$main_plots$scmet_mu_meth,
        bet_rds$main_plots$scmet_gamma_meth,
        wit_rds$main_plots$scmet_gamma_meth,
        bet_rds$main_plots$scmet_epsilon_meth,
        wit_rds$main_plots$scmet_epsilon_meth,
        legend = "none",
        ncol = 2,
        nrow = 3,
        heights = c(1, 1, 1)
      ),
      ggarrange(
        ggarrange(
          bet_rds$s_plots$secmem + xlab("# cells") + ylab("CPU time (s)"),
          bet_rds$n_plots$secmem + xlab("# features") + ylab("CPU time (s)"),
          wit_rds$s_plots$secmem + xlab("# cells") + ylab("CPU time (s)"),
          wit_rds$n_plots$secmem + xlab("# features") + ylab("CPU time (s)"),
          legend = "none",
          ncol = 4
        ),
        ggarrange(
          bet_rds$s_plots$nmi + xlab("# cells") + ylab("performance (metric NMI)"),
          bet_rds$n_plots$nmi + xlab("# features") + ylab("performance (metric NMI)"),
          wit_rds$s_plots$nmi + xlab("# cells") + ylab("performance (metric NMI)"),
          wit_rds$n_plots$nmi + xlab("# features") + ylab("performance (metric NMI)"),
          legend = "none",
          ncol = 4
        ),
        nrow = 2
      ),
      ncol = 2,
      widths = c(1, 2)
    ),
    nrow = 2,
    heights = c(1, 2.5)
  ),

  # shared legends stacked
  ggarrange(
    wit_rds$legends$bcVI,
    bet_rds$legends$wcVI,
    bet_rds$legends$nmi,
    bet_rds$legends$perf,
    ncol = 1,
    align = "v"
  ),
  ncol = 2,
  nrow = 1,
  widths = c(8, 1),
  align = "h"
)

ggsave(
  filename = snakemake@output[[1]],
  plot = arranged,
  width = 25,
  height = 15,
  units = "in"
)
