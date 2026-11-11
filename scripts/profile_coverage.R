# Coverage of profile-likelihood intervals for the variance ratio gamma in the
# replicate-trial experiment (weak, baseline, strong). Uses the same seeds,
# designs and data as run_replicate_trials() in R/experiments.R and saves
# results/profile_coverage.rds.
# From the repo root:   Rscript scripts/profile_coverage.R
# LVA_CORES (environment) sets the number of cores (default 2); N_REP overrides
# the replicate count for a quick dry run only.

source("R/settings.R")
for (f in setdiff(sort(list.files("R", pattern = "\\.R$", full.names = TRUE)),
                  "R/settings.R")) {
  source(f)
}

n_rep <- as.integer(Sys.getenv("N_REP", SETTINGS$n_rep))
n_cores <- as.integer(Sys.getenv("LVA_CORES", "2"))
max_cores <- max(1L, parallel::detectCores() - 1L)
if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L
n_cores <- min(n_cores, max_cores)
# mclapply forks, which Windows does not support: run serially there.
if (.Platform$OS.type == "windows") n_cores <- 1L

# lapply over the replicates, in parallel when possible. Every replicate sets
# its own seed, so the results do not depend on the route taken.
run_reps <- function(X, FUN) {
  if (n_cores > 1L) {
    out <- tryCatch(parallel::mclapply(X, FUN, mc.cores = n_cores),
                    error = function(e) {
                      message("mclapply failed (", conditionMessage(e),
                              "): falling back to serial")
                      NULL
                    })
    if (!is.null(out)) return(out)
  }
  lapply(X, FUN)
}

# Nominal 95% profile-likelihood interval: every gamma whose profile REML
# log-likelihood is within qchisq(0.95, 1) / 2 of the maximum.
DROP_95 <- qchisq(0.95, 1) / 2

# Interval for one dataset. The search runs over gamma = 0 plus the article's
# log-gamma grid; each end is refined with uniroot() between the last grid
# point inside and the first outside. If gamma = 0 is inside, the lower end is
# 0; if the top of the grid (gamma = 1e4) is inside, the upper end is Inf.
profile_interval <- function(setup, fit, gamma_true) {
  grid <- SETTINGS$log_gamma_grid
  cut <- fit$loglik - DROP_95
  ll_fun <- function(lg) reml_loglik(exp(lg), setup) - cut
  ll_grid <- reml_loglik(exp(grid), setup) - cut
  ll_zero <- reml_loglik(0, setup) - cut
  inside <- ll_grid >= 0
  # The peak itself may sit between grid points; it is always inside.
  lg_hat <- if (fit$boundary) -Inf else log(fit$gamma)

  # Number of separate pieces of the set on the grid (1 = an interval).
  pieces <- sum(diff(c(FALSE, c(ll_zero >= 0, inside))) == 1)

  if (ll_zero >= 0) {
    lo <- 0
  } else {
    i <- which(inside)[1]
    if (is.na(i)) i <- findInterval(lg_hat, grid) + 1   # only the peak is inside
    lo <- if (i == 1) {
      exp(uniroot(function(lg) ll_fun(lg), c(grid[1] - 30, grid[1]))$root)
    } else {
      upper <- if (inside[i]) grid[i] else lg_hat
      exp(uniroot(ll_fun, c(grid[i - 1], upper))$root)
    }
  }
  j <- tail(which(inside), 1)
  if (length(j) == 0) j <- findInterval(lg_hat, grid)
  hi <- if (j == length(grid)) Inf else {
    lower <- if (j >= 1 && inside[j]) grid[j] else lg_hat
    exp(uniroot(ll_fun, c(max(lower, if (j >= 1) grid[j] else -Inf), grid[j + 1]))$root)
  }
  ll_true <- reml_loglik(gamma_true, setup) - cut
  list(lo = lo, hi = hi, covered = ll_true >= 0,
       covered_interval = gamma_true >= lo && gamma_true <= hi,
       pieces = pieces)
}

one_rep <- function(scenario, r) {
  s <- SETTINGS[[scenario]]
  seed <- BASE_SEED + r
  design <- make_design(s$rep_counts, seed = seed)
  sim <- simulate_trial(design, h2 = s$h2, seed = seed)
  dm <- design_matrices(design)
  setup <- reml_setup(sim$y, dm$X, dm$Z)
  fit <- reml_fit(setup)
  gamma_true <- s$h2 / (1 - s$h2)
  iv <- profile_interval(setup, fit, gamma_true)
  data.frame(scenario = scenario, rep = as.integer(r),
             log_gamma_hat = log(fit$gamma), h2_hat = fit$h2,
             boundary = fit$boundary,
             gamma_lo = iv$lo, gamma_hi = iv$hi,
             h2_lo = iv$lo / (1 + iv$lo),
             h2_hi = if (is.finite(iv$hi)) iv$hi / (1 + iv$hi) else 1,
             reaches_zero = iv$lo == 0,
             covered = iv$covered, covered_interval = iv$covered_interval,
             pieces = iv$pieces, error = NA_character_)
}

t0 <- proc.time()[["elapsed"]]
out <- list()
for (sc in c("weak", "baseline", "strong")) {
  rows <- run_reps(seq_len(n_rep), function(r) {
    tryCatch(one_rep(sc, r), error = function(e) {
      data.frame(scenario = sc, rep = as.integer(r), log_gamma_hat = NA_real_,
                 h2_hat = NA_real_, boundary = NA, gamma_lo = NA_real_,
                 gamma_hi = NA_real_, h2_lo = NA_real_, h2_hi = NA_real_,
                 reaches_zero = NA, covered = NA, covered_interval = NA,
                 pieces = NA_integer_, error = conditionMessage(e))
    })
  })
  out[[sc]] <- do.call(rbind, rows)
  cat(sprintf("%-8s done at %6.1f s; failures: %d\n", sc,
              proc.time()[["elapsed"]] - t0, sum(!is.na(out[[sc]]$error))))
}
res <- do.call(rbind, out)
rownames(res) <- NULL

# Same estimates as the stored replicate experiment, replicate by replicate.
for (sc in names(out)) {
  stored <- readRDS(file.path("results", paste0("replicates_", sc, ".rds")))
  d <- out[[sc]]
  stored <- stored[match(d$rep, stored$rep), ]
  stopifnot(identical(d$h2_hat, stored$h2_hat),
            identical(d$log_gamma_hat, stored$log_gamma_hat),
            identical(d$boundary, stored$boundary))
}
cat("Estimates identical to results/replicates_*.rds\n")

attr(res, "n_rep") <- n_rep
attr(res, "drop") <- DROP_95
attr(res, "elapsed_s") <- proc.time()[["elapsed"]] - t0
saveRDS(res, file.path("results", "profile_coverage.rds"))

for (sc in names(out)) {
  d <- out[[sc]]
  cat(sprintf("%-8s covered %d / %d (%.1f%%); interval form %d; reaches 0: %d; split sets: %d\n",
              sc, sum(d$covered), nrow(d), 100 * mean(d$covered),
              sum(d$covered_interval), sum(d$reaches_zero), sum(d$pieces > 1)))
}
cat(sprintf("Elapsed: %.1f s on %d core(s)\n", attr(res, "elapsed_s"), n_cores))
