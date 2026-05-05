## Shared driver categorization for yamet reports.
##
## Classifies genomic annotations as "across-cell driven", "within-cell driven",
## "both", or "neither" based on how much adjH (across-cell heterogeneity) and
## adjS (within-cell heterogeneity) vary across biological groups.
##
## The logic: for each annotation, compute the SD of the group-level median
## adjH and adjS. If one SD is at least 1.5x the other, that entropy component
## dominates. If both SDs are below the 30th percentile of all annotations,
## the annotation is labelled "neither" (not varying enough to call).
## Otherwise it is "both".
##
## This never produces NAs: missing or zero SDs are mapped to "neither".
##
## Requires dplyr.

## Categorize annotations from a group-level summary data.frame.
##
## grp_df must contain columns: annotation, median_adjS, median_adjH, plus
## a grouping column named by `group_col` (e.g. "lineage_class", "cell_class",
## "location").
##
## Returns a data.frame with columns: annotation, adjH_sd, adjS_sd, driver.
categorize_drivers <- function(grp_df, group_col) {
  stopifnot(
    all(c("annotation", "median_adjS", "median_adjH", group_col) %in% names(grp_df))
  )

  grp_df <- grp_df[!is.na(grp_df[[group_col]]), , drop = FALSE]

  var_df <- grp_df %>%
    dplyr::group_by(annotation) %>%
    dplyr::summarise(
      adjH_sd = sd(median_adjH, na.rm = TRUE),
      adjS_sd = sd(median_adjS, na.rm = TRUE),
      .groups = "drop"
    )

  var_df$adjH_sd[!is.finite(var_df$adjH_sd)] <- 0
  var_df$adjS_sd[!is.finite(var_df$adjS_sd)] <- 0

  adjH_thr <- quantile(var_df$adjH_sd, 0.3)
  adjS_thr <- quantile(var_df$adjS_sd, 0.3)

  var_df$driver <- dplyr::case_when(
    var_df$adjH_sd < adjH_thr & var_df$adjS_sd < adjS_thr ~ "neither",
    var_df$adjH_sd >= var_df$adjS_sd * 1.5                 ~ "across-cell driven",
    var_df$adjS_sd >= var_df$adjH_sd * 1.5                 ~ "within-cell driven",
    TRUE                                                    ~ "both"
  )

  var_df
}

## Plot driver scatter: adjH_sd vs adjS_sd, colored/shaped by driver category.
## Returns a ggplot object. Requires ggplot2, ggrepel, and palettes.R loaded.
plot_driver_scatter <- function(driver_df, x_label = "SD of adjH across groups",
                                y_label = "SD of adjS across groups") {
  ggplot2::ggplot(driver_df,
                  ggplot2::aes(x = adjH_sd, y = adjS_sd,
                               color = driver, shape = driver,
                               label = annotation)) +
    ggplot2::geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = Inf,
                             box.padding = 0.5, point.padding = 0.2,
                             min.segment.length = 0, force = 2,
                             segment.size = 0.2, segment.alpha = 0.5) +
    ggplot2::scale_color_manual(values = driver_pal) +
    ggplot2::scale_shape_manual(values = driver_shapes) +
    ggplot2::labs(x = x_label, y = y_label, color = "driver", shape = "driver") +
    theme_ng() +
    ggplot2::theme(plot.margin = ggplot2::margin(3, 6, 3, 6, "mm"))
}
