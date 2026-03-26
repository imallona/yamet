library(ggplot2)
library(viridis)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

theme_ng <- function(base_size = 7, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid        = element_blank(),
      axis.line         = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks        = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = unit(1.5, "mm"),
      axis.text         = element_text(size = base_size, colour = "black"),
      axis.text.x       = element_text(size = base_size, colour = "black",
                                       angle = 45, hjust = 1, vjust = 1),
      axis.title        = element_text(size = base_size, colour = "black"),
      strip.background  = element_blank(),
      strip.text        = element_text(size = base_size, colour = "black"),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.key.size   = unit(3, "mm"),
      legend.text       = element_text(size = base_size - 1),
      legend.title      = element_text(size = base_size),
      plot.title        = element_text(size = base_size + 1, face = "plain",
                                       hjust = 0),
      plot.subtitle     = element_text(size = base_size, colour = "grey30"),
      plot.background   = element_blank(),
      complete          = TRUE
    )
}

## for scatter plots or UMAP axes where 45-degree rotation looks odd;
## silently drops labels that would overlap instead of rotating them
guide_x_nolap <- function() guide_axis(check.overlap = TRUE)

theme_ng_discrete <- function(base_size = 7, base_family = "Helvetica") {
  list(
    theme_ng(base_size = base_size, base_family = base_family),
    scale_color_manual(values = cbbPalette)
  )
}

theme_ng_continuous <- function(base_size = 7, base_family = "Helvetica") {
  list(
    theme_ng(base_size = base_size, base_family = base_family),
    scale_fill_viridis_c()
  )
}

## writes both PNG (600 dpi) and SVG; width/height in mm
save_ng <- function(plot, file, width_mm = 85, height_mm = 60) {
  ggsave(paste0(file, ".png"), plot,
         width = width_mm, height = height_mm, units = "mm", dpi = 600)
  ggsave(paste0(file, ".svg"), plot,
         width = width_mm, height = height_mm, units = "mm")
  invisible(plot)
}
