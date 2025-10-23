#' Differential entropy testing from a tidy data frame
#'
#' @param df       data.frame with columns: region (tile), sampen, meth, loc, patient (at minimum)
#' @param formula  model formula (e.g. sampen ~ meth + I(meth^2) + loc + patient), as string
#' @param param    BiocParallel parameter object
#' @param out_file path to save CSV results
#' @param top_n    number of top regions to return
#'
#' @return list with coefs_df, top_data (subset of df for top_n regions)
diff_entropy_test_df <- function(df,
                                 formula,
                                 param,
                                 out_file = "pt_vs_nc_coefs.csv",
                                 top_n = 2000,
                                 contrast = "locPT") {
    
   
  region_list <- split(df, df$region)
  
  rowwise_lm <- function(region_df) {
    fit <- try(lm(as.formula(formula), data = region_df), silent = TRUE)
    if (inherits(fit, "try-error")) return(rep(NA, 5))
    
    s <- summary(fit)$coefficients
    if (!(contrast %in% rownames(s))) return(rep(NA, 5))
    
    c(estimate  = s[contrast, "Estimate"],
      std_error = s[contrast, "Std. Error"],
      t_value   = s[contrast, "t value"],
      p_value   = s[contrast, "Pr(>|t|)"],
      df        = df.residual(fit))
  }
  
  coefs_list <- bplapply(region_list, rowwise_lm, BPPARAM = param)
  coefs_df <- as.data.frame(do.call(rbind, coefs_list))
  colnames(coefs_df) <- c("estimate", "std_error", "t_value", "p_value", "df")
  coefs_df$region <- names(region_list)
  coefs_df[] <- lapply(coefs_df, as.numeric)
  
  # Moderated t/p
  valid <- complete.cases(coefs_df)
  coefs_valid <- coefs_df[valid, ]
  
  s2 <- coefs_valid$std_error^2
  df_resid <- coefs_valid$df
  squeezed <- squeezeVar(var = s2, df = df_resid)
  
  moderated_t <- coefs_valid$estimate / sqrt(squeezed$var.post)
  moderated_p <- 2 * pt(-abs(moderated_t), df = squeezed$df.prior + df_resid)
  adj_p <- p.adjust(moderated_p, method = "BH")
  
  coefs_valid <- cbind(coefs_valid,
                       moderated_t = moderated_t,
                       moderated_p = moderated_p,
                       adj_p = adj_p)
  
  # merge back
  coefs_df <- merge(coefs_df,
                    coefs_valid[, c("region","moderated_t","moderated_p","adj_p")],
                    by = "region", all.x = TRUE)
  
  # output as CSV
  write.csv(coefs_df, file = out_file)
  
  # top hits
  sorted_idx <- order(coefs_df$adj_p, na.last = NA)
  top_idx <- head(sorted_idx, top_n)
  top_data <- df[df$region %in% coefs_df$region[top_idx], ]
  
  list(
    coefs_df = coefs_df,
    top_data = top_data
  )
}
