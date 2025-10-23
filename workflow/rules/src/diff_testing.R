#' Differential entropy testing with flexible formula and locations
#'
#' @param sub_sampens matrix of SampEn values (rows = regions, cols = cells)
#' @param sub_meths   matrix of methylation values (same dims as sub_sampens)
#' @param groups      data.frame with columns: subloc, patient (length = ncol)
#' @param formula     model formula (string or formula), e.g. "sampen ~ meth + I(meth^2) + loc + patient"
#' @param loc_levels  character vector of two location codes to compare, e.g. c("PT","NC")
#' @param ref_level   reference level for loc factor, e.g. "NC"
#' @param contrast    coefficient name to extract, e.g. "locPT"
#' @param param       BiocParallel parameter object
#' @param top_n       number of top regions to return (this is probably unnecessary)
#' @param out_file    file path to save coefs/ results (RDS)
#'
#' @return list with coefs_df, top_entropy, top_meth
diff_entropy_test <- function(sub_sampens, sub_meths, groups,
                              formula = "sampen ~ meth + I(meth^2) + loc + patient",
                              loc_levels = c("PT","NC"),
                              ref_level = "NC",
                              contrast = "locPT",
                              param,
                              top_n = 2000,
                              out_file = "diff_entropy_coefs.rds") {
  
  if (is.character(formula)) formula <- as.formula(formula)
  
  rowwise_lm <- function(i) {
    df <- data.frame(
      sampen  = sub_sampens[i, ],
      meth    = sub_meths[i, ],
      loc     = factor(substr(groups$subloc, 1, 2)),
      patient = groups$patient
    )
    df <- df[df$loc %in% loc_levels, ]
    df$loc <- relevel(factor(df$loc), ref = ref_level)
    
    fit <- try(lm(formula, data = df), silent = TRUE)
    if (inherits(fit, "try-error")) return(rep(NA, 5))
    
    s <- summary(fit)$coefficients
    if (!(contrast %in% rownames(s))) return(rep(NA, 5))
    
    c(
      estimate  = s[contrast, "Estimate"],
      std_error = s[contrast, "Std. Error"],
      t_value   = s[contrast, "t value"],
      p_value   = s[contrast, "Pr(>|t|)"],
      df        = df.residual(fit)
    )
  }
  
   coefs_list <- BiocParallel::bplapply(seq_len(nrow(sub_sampens)), rowwise_lm, BPPARAM = param)
  coefs_df <- as.data.frame(do.call(rbind, coefs_list))
  colnames(coefs_df) <- c("estimate","std_error","t_value","p_value","df")
  coefs_df[] <- lapply(coefs_df, as.numeric)
  
  # Moderated statistics
  valid <- complete.cases(coefs_df)
  coefs_valid <- coefs_df[valid, ]
  s2 <- coefs_valid$std_error^2
  df_resid <- coefs_valid$df
  squeezed <- limma::squeezeVar(var = s2, df = df_resid)
  
  moderated_t <- coefs_valid$estimate / sqrt(squeezed$var.post)
  moderated_p <- 2 * pt(-abs(moderated_t), df = squeezed$df.prior + df_resid)
  adj_p <- p.adjust(moderated_p, method = "BH")
  
  coefs_valid <- cbind(coefs_valid,
                       moderated_t = moderated_t,
                       moderated_p = moderated_p,
                       adj_p = adj_p)
  
  coefs_df$moderated_t <- NA
  coefs_df$moderated_p <- NA
  coefs_df$adj_p <- NA
  coefs_df[valid, c("moderated_t","moderated_p","adj_p")] <-
    coefs_valid[, c("moderated_t","moderated_p","adj_p")]
  
  coefs_df$region <- rownames(sub_sampens)
  
  saveRDS(coefs_df, file = out_file)
  
  ## top, even though they'd be filtered by significance
  sorted_idx <- order(coefs_df$adj_p, na.last = NA)
  top_idx <- head(sorted_idx, top_n)
  top_entropy <- sub_sampens[top_idx, , drop = FALSE]
  top_meth    <- sub_meths[top_idx, , drop = FALSE]
  
  list(
    coefs_df    = coefs_df,
    top_entropy = top_entropy,
    top_meth    = top_meth
  )
}
