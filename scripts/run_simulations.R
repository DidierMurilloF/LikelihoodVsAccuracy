# Run every experiment and save the results the article reads (results/*.rds).
# From the repo root:   Rscript scripts/run_simulations.R
# Replicate r always uses seed BASE_SEED + r, so the parallel run gives exactly
# the same results as a serial one (checked below before the main run).

source("R/settings.R")
for (f in setdiff(sort(list.files("R", pattern = "\\.R$", full.names = TRUE)),
                  "R/settings.R")) {
  source(f)
}
dir.create("results", showWarnings = FALSE)

# N_REP (environment) overrides the replicate count, for a quick dry run only.
n_rep <- as.integer(Sys.getenv("N_REP", SETTINGS$n_rep))
# LVA_CORES (environment) sets the number of cores (default 2), capped at
# detectCores() - 1 so the machine stays usable.
n_cores <- as.integer(Sys.getenv("LVA_CORES", "2"))
max_cores <- max(1L, parallel::detectCores() - 1L)
if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L
if (n_cores > max_cores) {
  message("LVA_CORES = ", n_cores, " is above detectCores() - 1; using ", max_cores)
  n_cores <- max_cores
}
cat(sprintf("Cores requested (LVA_CORES, default 2): %d\n", n_cores))

# Parallel sanity check: fall back to serial if mclapply misbehaves here.
if (n_cores > 1L) {
  n_cores <- tryCatch({
    same <- identical(run_replicate_trials("baseline", 16, cores = n_cores),
                      run_replicate_trials("baseline", 16))
    if (same) n_cores else {
      message("Parallel run differs from serial: falling back to serial")
      1L
    }
  }, error = function(e) {
    message("mclapply failed (", conditionMessage(e), "): falling back to serial")
    1L
  })
}
cat(sprintf("Running %d replicates per experiment on %d core(s)\n\n", n_rep, n_cores))

elapsed <- numeric()
failures <- integer()

# LVA_ONLY (environment, comma-separated experiment names, e.g.
# "two_analysts_spurious,two_analysts_real") re-runs only those experiments;
# the others are read from results/ (with their old elapsed time and failure
# count from summary.rds) and the summary is rebuilt from everything.
only <- Sys.getenv("LVA_ONLY", "")
only <- if (nzchar(only)) trimws(strsplit(only, ",")[[1]]) else NULL
old_summary <- if (!is.null(only) && file.exists("results/summary.rds")) {
  readRDS("results/summary.rds")
} else NULL

# Run one experiment, save it, and record its elapsed time and failures.
# `expr` is only evaluated when the experiment is selected.
run_and_save <- function(name, expr) {
  if (!is.null(only) && !name %in% only) {
    res <- readRDS(file.path("results", paste0(name, ".rds")))
    elapsed[[name]] <<- old_summary$elapsed_s[[name]]
    failures[[name]] <<- old_summary$failures[[name]]
    cat(sprintf("%-26s kept from results/\n", name))
    return(invisible(res))
  }
  t0 <- proc.time()[["elapsed"]]
  res <- expr
  saveRDS(res, file.path("results", paste0(name, ".rds")))
  secs <- proc.time()[["elapsed"]] - t0
  # Failed replicates: NA rows with an `error` message, or (curves) the
  # n_failed attribute.
  count_fail <- function(x) {
    if (is.data.frame(x) && "error" %in% names(x)) return(sum(!is.na(x$error)))
    if (is.data.frame(x)) return(as.integer(if (is.null(attr(x, "n_failed"))) 0 else attr(x, "n_failed")))
    sum(vapply(x, count_fail, numeric(1)))
  }
  n_fail <- count_fail(res)
  elapsed[[name]] <<- secs
  failures[[name]] <<- n_fail
  cat(sprintf("%-26s %7.1f s   failures: %d\n", name, secs, n_fail))
  invisible(res)
}

h2_kin <- SETTINGS$kinship$h2
names_kin <- paste0("h2_", h2_kin)

spurious <- run_and_save("two_analysts_spurious",
                         run_two_analysts("spurious_family", n_rep, cores = n_cores))
real <- run_and_save("two_analysts_real",
                     run_two_analysts("real_block", n_rep, cores = n_cores))
reps <- list()
for (sc in c("weak", "baseline", "strong")) {
  reps[[sc]] <- run_and_save(paste0("replicates_", sc),
                             run_replicate_trials(sc, n_rep, cores = n_cores))
}
run_and_save("curves_iid", run_curves("baseline"))
run_and_save("curves_kinship", setNames(
  lapply(h2_kin, function(h) run_curves("kinship", kinship = TRUE, h2 = h)), names_kin))
run_and_save("expected_curves_iid",
             run_expected_curves("baseline", n_rep = n_rep, cores = n_cores))
run_and_save("expected_curves_kinship", setNames(
  lapply(h2_kin, function(h) run_expected_curves("kinship", kinship = TRUE, h2 = h,
                                                 n_rep = n_rep, cores = n_cores)),
  names_kin))
conv <- run_and_save("convergence", run_convergence(n_rep, cores = n_cores))

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

# Two analysts. A's model nests B's, so loglik_A >= loglik_B. When A's extra
# variance is estimated as exactly 0 the two models coincide and the loglik
# difference is pure rounding (~1e-13, either sign), so "A higher" means
# loglik_A - loglik_B > LOGLIK_TOL. Accuracy comparisons use the replicates
# where both predictions are non-degenerate (cor defined); the degenerate
# counts are reported.
LOGLIK_TOL <- 1e-6
analysts_summary <- function(d) {
  d <- d[is.na(d$error), ]
  higher <- d$loglik_A - d$loglik_B > LOGLIK_TOL
  both_ok <- !d$degenerate_A & !d$degenerate_B
  worse <- d$cor_A < d$cor_B
  # 5% critical value of the 0.5 chi2_0 + 0.5 chi2_1 mixture.
  signif <- d$lrt > qchisq(0.90, 1)
  list(n = nrow(d),
       loglik_tol = LOGLIK_TOL,
       n_A_higher = sum(higher),
       p_A_higher = mean(higher),
       n_A_higher_raw = sum(d$loglik_A > d$loglik_B),
       n_extra_positive = sum(d$sigma2_extra_hat > 0),
       n_degenerate_A = sum(d$degenerate_A),
       n_degenerate_B = sum(d$degenerate_B),
       n_A_higher_both_ok = sum(higher & both_ok),
       p_A_worse_given_higher = mean(worse[higher & both_ok]),
       n_lrt_signif = sum(signif),
       n_lrt_signif_both_ok = sum(signif & both_ok),
       n_A_worse_given_signif = sum(worse[signif & both_ok]),
       mean_cor_A = mean(d$cor_A[both_ok]), mean_cor_B = mean(d$cor_B[both_ok]))
}

replicates_summary <- function(d) {
  ok <- is.na(d$error)
  ct <- cor.test(d$loglik_max, d$cor)          # complete cases (non-degenerate)
  list(n = sum(ok), n_boundary = sum(d$boundary[ok]),
       n_cor = sum(complete.cases(d$loglik_max, d$cor)),
       cor_loglik_acc = unname(ct$estimate),
       ci_low = ct$conf.int[1], ci_high = ct$conf.int[2])
}

res_summary <- list(
  n_rep = n_rep,
  cores = n_cores,
  analysts = list(spurious_family = analysts_summary(spurious),
                  real_block = analysts_summary(real)),
  replicates = lapply(reps, replicates_summary),
  convergence_maxit30_converged = with(conv[which(conv$maxit == 30), ],
                                       tapply(converged, start, mean)),
  failures = failures,
  elapsed_s = elapsed
)
saveRDS(res_summary, file.path("results", "summary.rds"))

cat("\n==== Summary ====\n")
cat("\nBoundary fits (sigma2_g_hat = 0) per scenario:\n")
for (sc in names(res_summary$replicates)) {
  s <- res_summary$replicates[[sc]]
  cat(sprintf("  %-8s %3d / %d (%.1f%%)\n", sc, s$n_boundary, s$n,
              100 * s$n_boundary / s$n))
}
cat("\nTwo analysts:\n")
for (tr in names(res_summary$analysts)) {
  s <- res_summary$analysts[[tr]]
  cat(sprintf(paste0("  %-15s P(loglik_A > loglik_B + %g) = %.3f (%d/%d; raw '>' %d; sigma2_extra > 0: %d)\n",
                     "  %-15s degenerate predictions: A %d, B %d\n",
                     "  %-15s when A higher, A less accurate (cor): %.1f%% (of %d)\n",
                     "  %-15s LRT 'significant' (> %.2f): %d, A less accurate in %d of %d\n",
                     "  %-15s mean cor (both non-degenerate): A %.4f, B %.4f\n"),
              tr, s$loglik_tol, s$p_A_higher, s$n_A_higher, s$n, s$n_A_higher_raw,
              s$n_extra_positive,
              "", s$n_degenerate_A, s$n_degenerate_B,
              "", 100 * s$p_A_worse_given_higher, s$n_A_higher_both_ok,
              "", qchisq(0.90, 1), s$n_lrt_signif, s$n_A_worse_given_signif,
              s$n_lrt_signif_both_ok,
              "", s$mean_cor_A, s$mean_cor_B))
}
cat("\ncor(loglik_max, cor) across replicate trials (95% CI):\n")
for (sc in names(res_summary$replicates)) {
  s <- res_summary$replicates[[sc]]
  cat(sprintf("  %-8s %.3f [%.3f, %.3f]  (n = %d non-degenerate)\n",
              sc, s$cor_loglik_acc, s$ci_low, s$ci_high, s$n_cor))
}
cat("\nConvergence, share converged at maxit = 30:\n")
print(res_summary$convergence_maxit30_converged)
cat("\nFailed replicates per experiment:\n")
print(unlist(failures))
cat(sprintf("\nTotal elapsed: %.1f s\n", sum(unlist(elapsed))))
