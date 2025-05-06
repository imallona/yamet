#' Helper function to generate a table of normalised mutual
#' information scores from true and computed heterogeneities
#'
#' @author Atreya Choudhury
#' @date 2025-04-23

library(infotheo)
library(knitr)

mi_table_gen <- function(jnt, scmet_fit, truth, table = F) {
  te <- entropy(truth)
  mi_sampen <- mutinformation(discretize(jnt$sampen_avg), truth) /
    sqrt(entropy(discretize(jnt$sampen_avg)) * te)
  mi_shannon <- mutinformation(discretize(jnt$shannon), truth) /
    sqrt(entropy(discretize(jnt$shannon)) * te)
  mi_mu <- mutinformation(discretize(scmet_fit$posterior$mu_avg), truth) /
    sqrt(entropy(discretize(scmet_fit$posterior$mu_avg)) * te)
  mi_gamma <- mutinformation(discretize(scmet_fit$posterior$gamma_avg), truth) /
    sqrt(entropy(discretize(scmet_fit$posterior$gamma_avg)) * te)
  mi_epsilon <- mutinformation(discretize(scmet_fit$posterior$epsilon_avg), truth) /
    sqrt(entropy(discretize(scmet_fit$posterior$epsilon_avg)) * te)

  mi_results <- data.frame(
    "Metric" = c(
      "Sample Entropy", "Shannon Entropy",
      "scMET (mu)", "scMET (gamma)", "scMET (epsilon)"
    ),
    "NMI" = c(
      mi_sampen, mi_shannon, mi_mu, mi_gamma, mi_epsilon
    ),
    check.names = FALSE
  )

  if (table) {
    kable(mi_results,
      caption = "Normalized Mutual Information with True Heterogeneity Index",
      col.names = c("Metric", "Normalised Mutual Information (NMI)"),
      digits = 4,
      align = c("l", "r")
    )
  } else {
    return(mi_results)
  }
}
