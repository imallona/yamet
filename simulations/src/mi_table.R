#' Helper function to generate a table of normalised mutual
#' information scores from true and computed heterogeneities
#'
#' @author Atreya Choudhury
#' @date 2025-04-23

library(infotheo)
library(ggpubr)

mi_table_gen <- function(jnt, index_name, table = F) {
  index <- jnt[[index_name]]
  te <- entropy(index)

  calc_nmi <- function(x, index, te) {
    disc_x <- discretize(x)
    mi <- mutinformation(disc_x, index)
    hx <- entropy(disc_x)
    mi / sqrt(hx * te)
  }

  nmi_values <- c(
    "Sample Entropy" = calc_nmi(jnt$sampen_avg, index, te),
    "Shannon Entropy" = calc_nmi(jnt$shannon, index, te),
    "Average Methylation" = calc_nmi(jnt$avg_meth, index, te),
    "scMET (mu)" = calc_nmi(jnt$scmet_mu, index, te),
    "scMET (gamma)" = calc_nmi(jnt$scmet_gamma, index, te),
    "scMET (epsilon)" = calc_nmi(jnt$scmet_epsilon, index, te)
  )

  mi_results <- data.frame(
    Metric = names(nmi_values),
    NMI = unname(nmi_values),
    stringsAsFactors = FALSE
  )

  if (table) {
    tbl <- ggtexttable(
      mi_results,
      rows = NULL,
      cols = c("Metric", "Normalized Mutual Information (NMI)"),
      theme = ttheme(
        colnames.style = colnames_style(
          fill = "white"
        ),
        tbody.style = tbody_style(
          fill = "white",
          hjust = as.vector(matrix(c(0, 1), ncol = 2, nrow = nrow(mi_results), byrow = TRUE)),
          x = as.vector(matrix(c(.1, .9), ncol = 2, nrow = nrow(mi_results), byrow = TRUE))
        )
      )
    )
    tbl <- tab_add_hline(tbl, at.row = 1:2, row.side = "top", linewidth = 2)
    return(tbl)
  } else {
    return(mi_results)
  }
}
