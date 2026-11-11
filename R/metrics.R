# Accuracy and scale metrics: predicted genetic values scored against the
# simulated true values.

# TRUE when g_hat carries no ranking information (all values equal).
is_degenerate <- function(g_hat) {
  isTRUE(sd(g_hat) == 0) || all(g_hat == g_hat[1])
}

# Accuracy of g_hat against g_true, relative to the oracle predictions.
#   cor     Pearson correlation of g_hat with g_true
#   topk    overlap of the top-k sets of g_hat and g_true, as a fraction of k
#           (k = round(top_frac * n), at least 1; ties broken by position)
#   d_cor   cor(g_hat) - cor(g_oracle)
#   d_mse   MSE(g_hat) - MSE(g_oracle), MSE = mean((g - g_true)^2)
# When g_hat is constant (e.g. the gamma = 0 boundary), cor, topk and d_cor
# are NA_real_ and degenerate is TRUE; d_mse is still computed.
accuracy_metrics <- function(g_hat, g_true, g_oracle, top_frac = SETTINGS$top_frac) {
  mse <- function(x) mean((x - g_true)^2)
  d_mse <- mse(g_hat) - mse(g_oracle)
  degenerate <- is_degenerate(g_hat)
  if (degenerate) {
    return(list(cor = NA_real_, topk = NA_real_, d_cor = NA_real_,
                d_mse = d_mse, degenerate = TRUE))
  }
  k <- max(1L, as.integer(round(top_frac * length(g_true))))
  top_set <- function(x) order(-x)[seq_len(k)]
  r <- cor(g_hat, g_true)
  r_oracle <- if (is_degenerate(g_oracle)) NA_real_ else cor(g_oracle, g_true)
  list(cor = r,
       topk = length(intersect(top_set(g_hat), top_set(g_true))) / k,
       d_cor = r - r_oracle,
       d_mse = d_mse,
       degenerate = FALSE)
}

# Scale diagnostics for g_hat.
#   gain_slope   slope of lm(g_true ~ g_hat) (1 = well calibrated; > 1 means
#                the predictions are shrunk too much)
#   reliability  mean(1 - pev_scaled * sigma2_e / sigma2_g), the model's own
#                claim about its accuracy squared
# Both are NA_real_ when g_hat is constant; reliability is also NA_real_ when
# pev_scaled is NA or sigma2_g is 0.
scale_metrics <- function(g_hat, g_true, pev_scaled, sigma2_g, sigma2_e) {
  if (is_degenerate(g_hat)) {
    return(list(gain_slope = NA_real_, reliability = NA_real_))
  }
  gain_slope <- unname(coef(lm(g_true ~ g_hat))[2])
  reliability <- if (anyNA(pev_scaled) || !(sigma2_g > 0)) NA_real_ else
    mean(1 - pev_scaled * sigma2_e / sigma2_g)
  list(gain_slope = gain_slope, reliability = reliability)
}
