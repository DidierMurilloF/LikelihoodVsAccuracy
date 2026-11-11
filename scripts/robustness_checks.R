# Two robustness checks for the Myth 2 claim that, in a simple genotype model,
# rankings barely depend on the variance estimate. Saves
# results/robustness_checks.rds.
#
# (a) Uneven replication: 10 genotypes with 8 plots and 190 with 1 plot
#     (h2 = 0.3, otherwise the baseline set-up), with the baseline design on
#     the same seeds as a control.
# (b) Kinship within and between families: the article's genomic G at
#     h2 = 0.3 and 0.1, with accuracy computed over all genotypes, within
#     families (mean of the 10 within-family correlations) and between
#     families (correlation of the family means).
#
# In each trial the BLUPs are computed at the true variance ratio, at ratios
# 10 times too small and too large (the article's Myth 2 cost), and at the REML
# estimate. Replicate r uses seed BASE_SEED + r, as every experiment in
# R/experiments.R does, so the kinship trials are the first N_REP trials of the
# article's kinship experiment.
#
# From the repo root:   Rscript scripts/robustness_checks.R
# LVA_CORES (environment) sets the number of cores (default 2); N_REP sets the
# number of trials per setting (default 150).

source("R/settings.R")
for (f in setdiff(sort(list.files("R", pattern = "\\.R$", full.names = TRUE)),
                  "R/settings.R")) {
  source(f)
}

n_rep <- as.integer(Sys.getenv("N_REP", "150"))
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

WRONG <- 10                                     # "ten times too small / too large"
UNEVEN_REPS <- c(rep(8L, 10), rep(1L, 190))     # n = 270 plots
POINTS <- c("lo", "true", "hi", "reml")

# The settings checked. `family` labels genotypes for the within/between
# split (kinship only); `high_rep` marks the heavily replicated genotypes.
G_kin <- .kinship_G()
m_kin <- nrow(G_kin)
fam_kin <- ceiling(seq_len(m_kin) * SETTINGS$kinship$n_family / m_kin)
checks <- list(
  uneven = list(rep_counts = UNEVEN_REPS, h2 = 0.3, G = NULL, family = NULL),
  baseline = list(rep_counts = SETTINGS$baseline$rep_counts,
                  h2 = SETTINGS$baseline$h2, G = NULL, family = NULL),
  kinship_0.3 = list(rep_counts = SETTINGS$kinship$rep_counts,
                     h2 = SETTINGS$kinship$h2[1], G = G_kin, family = fam_kin),
  kinship_0.1 = list(rep_counts = SETTINGS$kinship$rep_counts,
                     h2 = SETTINGS$kinship$h2[2], G = G_kin, family = fam_kin)
)

# Within-family accuracy (mean of the per-family correlations) and
# between-family accuracy (correlation of family means).
family_cors <- function(gh, g, fam) {
  if (is.null(fam) || is_degenerate(gh)) return(c(within = NA_real_, between = NA_real_))
  within <- vapply(split(seq_along(g), fam), function(i) cor(gh[i], g[i]), numeric(1))
  c(within = mean(within),
    between = cor(tapply(gh, fam, mean), tapply(g, fam, mean)))
}

one_rep <- function(ck, r) {
  seed <- BASE_SEED + r
  design <- make_design(ck$rep_counts, seed = seed)
  sim <- simulate_trial(design, h2 = ck$h2, G = ck$G, seed = seed)
  dm <- design_matrices(design)
  fit <- reml_fit(reml_setup(sim$y, dm$X, dm$Z, ck$G))
  g0 <- ck$h2 / (1 - ck$h2)
  gammas <- c(g0 / WRONG, g0, g0 * WRONG)
  g_hat <- .grid_predictions(sim$y, dm$X, dm$Z, ck$G, gammas)
  # At the boundary every BLUP is 0: no ranking (degenerate).
  g_reml <- if (fit$boundary) rep(0, length(sim$g)) else
    drop(.grid_predictions(sim$y, dm$X, dm$Z, ck$G, fit$gamma))
  g_hat <- cbind(g_hat, g_reml)
  high <- ck$rep_counts > 1 & ck$rep_counts == max(ck$rep_counts)
  k <- max(1L, as.integer(round(SETTINGS$top_frac * length(sim$g))))
  rows <- lapply(seq_along(POINTS), function(j) {
    gh <- g_hat[, j]
    acc <- accuracy_metrics(gh, sim$g, g_hat[, 2])
    fc <- family_cors(gh, sim$g, ck$family)
    top <- if (acc$degenerate) integer() else order(-gh)[seq_len(k)]
    data.frame(rep = as.integer(r), point = POINTS[j],
               gamma = if (j == 4) fit$gamma else gammas[j],
               boundary = fit$boundary, cor = acc$cor, topk = acc$topk,
               cor_within = fc[["within"]], cor_between = fc[["between"]],
               n_high_in_top = if (acc$degenerate) NA_integer_ else sum(high[top]),
               degenerate = acc$degenerate, error = NA_character_)
  })
  do.call(rbind, rows)
}

t0 <- proc.time()[["elapsed"]]
out <- list()
for (nm in names(checks)) {
  ck <- checks[[nm]]
  rows <- run_reps(seq_len(n_rep), function(r) {
    tryCatch(one_rep(ck, r), error = function(e) {
      data.frame(rep = as.integer(r), point = NA_character_, gamma = NA_real_,
                 boundary = NA, cor = NA_real_, topk = NA_real_,
                 cor_within = NA_real_, cor_between = NA_real_,
                 n_high_in_top = NA_integer_, degenerate = NA,
                 error = conditionMessage(e))
    })
  })
  d <- do.call(rbind, rows)
  d <- cbind(check = nm, h2 = ck$h2, d)
  out[[nm]] <- d
  cat(sprintf("%-12s done at %6.1f s; failed trials: %d\n", nm,
              proc.time()[["elapsed"]] - t0, length(unique(d$rep[!is.na(d$error)]))))
}
res <- do.call(rbind, out)
rownames(res) <- NULL

# Relationships in G within and between families (off-diagonal elements).
same_fam <- outer(fam_kin, fam_kin, "==")
off <- row(G_kin) != col(G_kin)
attr(res, "G_within") <- mean(G_kin[same_fam & off])
attr(res, "G_between") <- mean(G_kin[!same_fam])
attr(res, "uneven_rep_counts") <- UNEVEN_REPS
attr(res, "wrong") <- WRONG
attr(res, "n_rep") <- n_rep
attr(res, "cores") <- n_cores
attr(res, "elapsed_s") <- proc.time()[["elapsed"]] - t0
saveRDS(res, file.path("results", "robustness_checks.rds"))

# Summary: mean of each metric at each point (non-degenerate trials).
ok <- is.na(res$error) & !res$degenerate
for (nm in names(checks)) {
  d <- res[ok & res$check == nm, ]
  cat(sprintf("\n%s (h2 = %.1f), %d trials, %d boundary fits\n", nm, checks[[nm]]$h2,
              length(unique(d$rep)), length(unique(res$rep[res$check == nm & res$boundary %in% TRUE]))))
  agg <- sapply(POINTS, function(p) {
    x <- d[d$point == p, ]
    c(cor = mean(x$cor), topk = mean(x$topk), within = mean(x$cor_within),
      between = mean(x$cor_between), high_top = mean(x$n_high_in_top))
  })
  print(round(agg, 4))
}
cat(sprintf("\nG: mean within-family relationship %.3f, between %.3f\n",
            attr(res, "G_within"), attr(res, "G_between")))
cat(sprintf("Elapsed: %.1f s on %d core(s)\n", attr(res, "elapsed_s"), n_cores))
