#!/usr/bin/env Rscript
## Figure 1, Panel A. Schematic of the single-cell DNA methylation data layout
## that yamet consumes.
##
## A toy 8-cell by 16-CpG matrix with:
##   - top track:  the genomic feature each CpG falls into
##                 (promoter, gene body, enhancer)
##   - left track: cell cluster membership (three clusters)
##   - matrix:     methylation status per cell and CpG, drawn as filled circles
##                 (methylated), open circles (unmethylated), or grey circles
##                 (missing / uncovered).
##
## Output:
##   <out_dir>/yamet_fig1_panel_a.svg   vector, drop into the Inkscape master
##   <out_dir>/yamet_fig1_panel_a.pdf   cairo_pdf, embed directly in LaTeX
##
## Usage:
##   Rscript fig1_panel_a.R [out_dir]
##
## Defaults to ./figures/ relative to the current working directory when no
## output directory is passed.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else file.path(getwd(), "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(2026)

n_cells <- 8L
n_cpgs  <- 16L

cluster <- factor(rep(c("Cluster 1", "Cluster 2", "Cluster 3"), c(3, 3, 2)),
                  levels = c("Cluster 1", "Cluster 2", "Cluster 3"))

feature <- factor(c(rep("Promoter", 5),
                    rep("Gene body", 8),
                    rep("Enhancer", 3)),
                  levels = c("Promoter", "Gene body", "Enhancer"))

make_row <- function(cluster_id) {
  x <- integer(n_cpgs)
  if (cluster_id == 1L) {
    x[1:5]   <- rbinom(5, 1, 0.85)
    x[6:13]  <- rbinom(8, 1, 0.5)
    x[14:16] <- rbinom(3, 1, 0.2)
  } else if (cluster_id == 2L) {
    x <- rbinom(n_cpgs, 1, 0.5)
  } else {
    x[1:5]   <- rbinom(5, 1, 0.2)
    x[6:13]  <- rbinom(8, 1, 0.5)
    x[14:16] <- rbinom(3, 1, 0.85)
  }
  x
}

M <- t(vapply(as.integer(cluster), make_row, integer(n_cpgs)))
miss <- matrix(rbinom(n_cells * n_cpgs, 1, 0.1) == 1L, nrow = n_cells)
M[miss] <- NA_integer_

df <- expand.grid(cell = seq_len(n_cells), cpg = seq_len(n_cpgs))
df$value <- as.vector(M)
df$state <- factor(
  ifelse(is.na(df$value), "missing",
         ifelse(df$value == 1L, "methylated", "unmethylated")),
  levels = c("methylated", "unmethylated", "missing")
)

clu_df  <- data.frame(cell = seq_len(n_cells), cluster = cluster)
feat_df <- data.frame(cpg  = seq_len(n_cpgs),  feature = feature)

state_fill <- c(
  "methylated"   = "#000000",
  "unmethylated" = "#FFFFFF",
  "missing"      = "#CCCCCC"
)
state_stroke <- c(
  "methylated"   = "#000000",
  "unmethylated" = "#000000",
  "missing"      = "#888888"
)
cluster_pal <- c(
  "Cluster 1" = "#1B9E77",
  "Cluster 2" = "#7570B3",
  "Cluster 3" = "#D95F02"
)
feature_pal <- c(
  "Promoter"  = "#377EB8",
  "Gene body" = "#4DAF4A",
  "Enhancer"  = "#984EA3"
)

base_theme <- theme_void(base_size = 7, base_family = "Helvetica") +
  theme(
    plot.margin   = unit(c(0, 0, 0, 0), "mm"),
    legend.title  = element_text(size = 7, colour = "black"),
    legend.text   = element_text(size = 6, colour = "black"),
    legend.key.size = unit(2.5, "mm")
  )

p_feat <- ggplot(feat_df, aes(x = cpg, y = 1, fill = feature)) +
  geom_tile(width = 0.96, height = 0.6) +
  scale_fill_manual(values = feature_pal, name = "Feature") +
  scale_x_continuous(limits = c(0.5, n_cpgs + 4.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  base_theme

p_clu <- ggplot(clu_df, aes(x = 1, y = cell, fill = cluster)) +
  geom_tile(width = 0.6, height = 0.96) +
  scale_fill_manual(values = cluster_pal, name = "Cluster") +
  scale_x_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  scale_y_reverse(limits = c(n_cells + 0.5, 0.5), expand = c(0, 0)) +
  base_theme

## Boxes showing the unit of computation:
##   adjS is per cell x per feature (one row x one feature span)
##   adjH is per cluster x per feature (rows of one cluster x one feature span)
## Both boxes pick the gene-body feature span (CpGs 6-13) so they can be
## compared directly: adjS uses cell 5 (cluster 2), adjH uses cluster 1.
adjS_box <- data.frame(xmin = 5.5, xmax = 13.5, ymin = 4.5, ymax = 5.5)
adjH_box <- data.frame(xmin = 5.5, xmax = 13.5, ymin = 0.5, ymax = 3.5)

p_main <- ggplot(df, aes(x = cpg, y = cell)) +
  geom_point(aes(fill = state, colour = state),
             shape = 21, size = 2.6, stroke = 0.4) +
  geom_rect(data = adjS_box, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, colour = "#222222", linewidth = 0.55,
            linetype = "dashed") +
  geom_rect(data = adjH_box, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, colour = "#222222", linewidth = 0.55,
            linetype = "solid") +
  annotate("text", x = 13.7, y = 5, hjust = 0, vjust = 0.5,
           label = "italic(S)~per~cell~x~feature",
           parse = TRUE, size = 2.3, colour = "#222222") +
  annotate("text", x = 13.7, y = 2, hjust = 0, vjust = 0.5,
           label = "italic(H)~per~cluster~x~feature",
           parse = TRUE, size = 2.3, colour = "#222222") +
  scale_fill_manual(values = state_fill,    name = "Methylation") +
  scale_colour_manual(values = state_stroke, guide = "none") +
  scale_x_continuous(limits = c(0.5, n_cpgs + 4.5), expand = c(0, 0),
                     breaks = NULL) +
  scale_y_reverse(limits = c(n_cells + 0.5, 0.5), expand = c(0, 0),
                  breaks = NULL) +
  base_theme +
  labs(x = "Cytosines along the genome", y = "Cells") +
  theme(
    axis.title.x = element_text(size = 7, margin = margin(t = 1)),
    axis.title.y = element_text(size = 7, angle = 90, margin = margin(r = 1))
  )

spacer <- plot_spacer()

panel <- (spacer + p_feat) / (p_clu + p_main) +
  plot_layout(widths  = c(1, 12),
              heights = c(1, 8),
              guides  = "collect") &
  theme(legend.position = "right",
        legend.box      = "vertical",
        legend.spacing.y = unit(2, "mm"))

svg_path <- file.path(out_dir, "yamet_fig1_panel_a.svg")
pdf_path <- file.path(out_dir, "yamet_fig1_panel_a.pdf")

ggsave(svg_path, panel, width = 115, height = 55, units = "mm",
       device = svglite)
ggsave(pdf_path, panel, width = 115, height = 55, units = "mm",
       device = cairo_pdf)

message("Wrote: ", svg_path)
message("Wrote: ", pdf_path)
