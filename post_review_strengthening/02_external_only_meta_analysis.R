# TNBC-Triclosan post-rejection strengthening
# External-only meta-analysis: GSE31519 + GSE58812 only
# Purpose: test whether the central EMT association persists after removing the
# TCGA derivation/internal-association cohort. This is a robustness analysis,
# not a search for a new positive result.

rm(list = ls(all.names = TRUE))
options(stringsAsFactors = FALSE, warn = 1)

INPUT <- "01_external_only_input_effects.csv"
OUT <- "03_external_only_meta_results.csv"

stopf <- function(...) stop(sprintf(...), call. = FALSE)
assert <- function(x, ...) if (!isTRUE(x)) stopf(...)
clip_r <- function(x) pmax(pmin(as.numeric(x), 0.999999), -0.999999)

x <- read.csv(INPUT, check.names = FALSE)
required <- c("cohort","effect_set","sample_number","rho","fisher_z",
              "empirical_fisher_z_variance","conventional_fisher_z_variance")
assert(all(required %in% names(x)), "Input columns are incomplete.")
assert(nrow(x) == 4L, "Expected four effect rows, observed %d.", nrow(x))
assert(setequal(unique(x$cohort), c("GSE31519","GSE58812")),
       "Only the two independent external cohorts are allowed.")
assert(setequal(unique(x$effect_set), c("FULL_EMT","NO_OVERLAP_EMT")),
       "Effect-set labels are incorrect.")
assert(all(x$sample_number %in% c(579L,107L)), "Sample-number lock failed.")
assert(all(abs(x$fisher_z - atanh(clip_r(x$rho))) < 1e-12),
       "Fisher-z values do not reproduce rho.")

meta_one <- function(d, variance_col, variance_method) {
  z <- as.numeric(d$fisher_z)
  v <- as.numeric(d[[variance_col]])
  assert(length(z) == 2L && all(is.finite(v) & v > 0),
         "Invalid external-only inputs.")
  k <- length(z)
  w <- 1 / v
  z_fixed <- sum(w * z) / sum(w)
  se_fixed <- sqrt(1 / sum(w))
  Q <- sum(w * (z - z_fixed)^2)
  Q_df <- k - 1L
  Q_p <- pchisq(Q, df = Q_df, lower.tail = FALSE)
  I2 <- if (Q > 0) max(0, (Q - Q_df) / Q) * 100 else 0
  C <- sum(w) - sum(w^2) / sum(w)
  tau2 <- if (C > 0) max(0, (Q - Q_df) / C) else 0
  wr <- 1 / (v + tau2)
  z_random <- sum(wr * z) / sum(wr)
  se_random <- sqrt(1 / sum(wr))
  scale <- max(1, sum(wr * (z - z_random)^2) / (k - 1L))
  se_mkh <- sqrt(scale / sum(wr))

  one_row <- function(method, est, se, dist = "NORMAL") {
    if (dist == "NORMAL") {
      crit <- qnorm(0.975)
      p <- 2 * pnorm(abs(est / se), lower.tail = FALSE)
      df <- NA_real_
    } else {
      df <- k - 1L
      crit <- qt(0.975, df = df)
      p <- 2 * pt(abs(est / se), df = df, lower.tail = FALSE)
    }
    lo <- est - crit * se
    hi <- est + crit * se
    data.frame(
      effect_set = unique(d$effect_set),
      variance_method = variance_method,
      pooling_method = method,
      cohort_number = k,
      total_sample_number = sum(d$sample_number),
      pooled_fisher_z = est,
      standard_error = se,
      inference_distribution = dist,
      inference_df = df,
      pooled_rho = tanh(est),
      confidence_interval_lower_rho = tanh(lo),
      confidence_interval_upper_rho = tanh(hi),
      two_sided_p = p,
      Cochran_Q = Q,
      Q_df = Q_df,
      Q_p = Q_p,
      I2_percent = I2,
      tau2 = tau2,
      modified_KH_residual_scale = if (method == "RANDOM_DL_MODIFIED_KH") scale else NA_real_,
      all_external_directions_positive = all(d$rho > 0),
      stringsAsFactors = FALSE
    )
  }

  rbind(
    one_row("FIXED_INVERSE_VARIANCE_NORMAL", z_fixed, se_fixed, "NORMAL"),
    one_row("RANDOM_DL_NORMAL", z_random, se_random, "NORMAL"),
    one_row("RANDOM_DL_MODIFIED_KH", z_random, se_mkh, "T")
  )
}

out <- do.call(rbind, lapply(split(x, x$effect_set), function(d) {
  rbind(
    meta_one(d, "empirical_fisher_z_variance", "EMPIRICAL_BOOTSTRAP_FISHER_Z"),
    meta_one(d, "conventional_fisher_z_variance", "CONVENTIONAL_1_OVER_N_MINUS_3")
  )
}))

write.csv(out, OUT, row.names = FALSE, na = "")
cat("External-only meta-analysis complete.\n")
print(out)
