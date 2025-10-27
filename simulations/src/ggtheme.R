library(viridis)

# Taken from: https://ghurault.github.io/HuraultMisc/reference/cbbPalette.html
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

theme_ng <- function(base_size = 12, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      # Remove grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      # Keep full axes
      axis.line = element_line(color = "black", size = 0.5),

      # Ticks and text
      axis.ticks = element_line(color = "black", size = 0.5),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold"),

      # Legend
      legend.key = element_blank(),
      legend.text = element_text(size = base_size - 1),
      legend.title = element_text(face = "bold"),

      # Strip text for facets
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),

      # Plot title
      plot.title = element_text(face = "bold", hjust = 0.5),

      # Square aspect ratio
      aspect.ratio = 1
    )
}

theme_ng_discrete <- function(base_size = 12, base_family = "sans") {
  list(
    theme_ng(base_size = base_size, base_family = base_family),
    scale_color_manual(values=cbbPalette)
  )
}

theme_ng_continuous <- function(base_size = 12, base_family = "sans") {
  list(
    theme_ng(base_size = base_size, base_family = base_family),
    scale_fill_viridis_d()
  )
}
