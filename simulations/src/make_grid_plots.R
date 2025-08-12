suppressPackageStartupMessages({
  library(ggplot2)
})

source("src/mi_table.R")
source("src/ggtheme.R")

get_combs <- function(GRID) {
  combinations <- expand.grid(GRID)
  tuples <- lapply(seq_len(nrow(combinations)), function(i) {
    as.list(combinations[i, ])
  })
  return(tuples)
}

get_mi <- function(data_dir, S, N, f, x) {
  intervals <- read.table(
    paste0(data_dir, paste("/intervals", S, N, f, "tsv", sep = ".")),
    header = T
  )
  intervals <- intervals[x]
  sampens <- read.table(
    paste0(data_dir, paste("/yamet", S, N, f, "det.out", sep = ".")),
    header = T
  )
  jnt <- cbind(sampens, intervals)
  rownames(jnt) <- paste(jnt$start, jnt$end, sep = ":")
  output_cols <- grep("^output", colnames(jnt), value = TRUE)

  scmet_fit <- readRDS(paste0(data_dir, paste("/scmet", S, N, f, "rds", sep = ".")))
  tmp <- colMeans(scmet_fit$posterior$mu)
  jnt$scmet_mu <- tmp[match(rownames(jnt), names(tmp))]
  tmp <- colMeans(scmet_fit$posterior$gamma)
  jnt$scmet_gamma <- tmp[match(rownames(jnt), names(tmp))]
  tmp <- colMeans(scmet_fit$posterior$epsilon)
  jnt$scmet_epsilon <- tmp[match(rownames(jnt), names(tmp))]

  jnt$sampen_avg <- rowMeans(jnt[, output_cols])

  mi_table <- mi_table_gen(jnt, x)
  return(mi_table)
}

make_grid_plots <- function(GRID, grid_var, x, data_dir) {
  bm_data <- do.call(
    rbind,
    lapply(get_combs(GRID), function(row) {
      S <- row$S
      N <- row$N
      f <- row$f
      do.call(rbind, lapply(c("yamet", "scmet"), function(mtd) {
        transform(read.table(paste0(
          data_dir, "/", paste(mtd, S, N, f, "benchmark", sep = ".")
        ), header = T), mtd = mtd, S = S, N = N, f = f)
      }))
    })
  )

  mi_data <- do.call(
    rbind,
    lapply(get_combs(GRID), function(row) {
      S <- row$S
      N <- row$N
      f <- row$f
      transform(get_mi(data_dir, S, N, f, x), S = S, N = N, f = f)
    })
  )

  plots <- list()

  plots$nmi <- ggplot(mi_data, aes(x = .data[[grid_var]], y = NMI, color = Metric, group = Metric)) +
    # scale_y_log10() +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    # scale_color_manual(values = c("yamet" = "#E69F00", "scmet" = "#56B4E9")) +
    labs(
      ## title = "NMI Comparison",
      x = paste(grid_var), y = "NMI",
      color = "Metric"
    ) +
    theme_ng() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )

  plots$seconds <- ggplot(bm_data, aes(x = .data[[grid_var]], y = s, color = mtd, group = mtd)) +
    scale_y_log10() +
    geom_point(size = 3, alpha = 0.7) +
    stat_summary(fun = mean, geom = "line", linewidth = 1.2, aes(group = mtd)) +
    scale_color_manual(values = c("yamet" = "#E69F00", "scmet" = "#56B4E9")) +
    labs(
      x = paste(grid_var), y = "seconds",
      color = "Method"
    ) +
    theme_ng() +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )

  plots$mem <- ggplot(bm_data, aes(x = .data[[grid_var]], y = max_rss, color = mtd, group = mtd)) +
    scale_y_log10() +
    geom_point(size = 3) +
    stat_summary(fun = mean, geom = "line", linewidth = 1.2, aes(group = mtd)) +
    scale_color_manual(values = c("yamet" = "#E69F00", "scmet" = "#56B4E9")) +
    labs(
      x = paste(grid_var), y = "MB",
      color = "Method"
    ) +
    theme_ng() +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )

  return(plots)
}
