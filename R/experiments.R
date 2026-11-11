# Experiments behind every number in the article. Each run_* function returns
# a plain data.frame; replicate r always uses seed BASE_SEED + r (for the field
# layout and for the simulated data), so serial and parallel runs give
# identical results.
#
# Robustness: every replicate runs inside tryCatch(). A failed replicate is
# kept as a row of NA with its message in the `error` column (NA when the
# replicate succeeded); rows are never dropped.

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

# Run one_rep(r) for every r in reps (in parallel when cores > 1) and stack
# the rows. `proto` is a zero-row data.frame with the columns of a successful
# row (without `error`); it gives a failed replicate its NA row.
.run_reps <- function(reps, one_rep, proto, cores = 1L) {
  na_row <- function(r, msg) {
    row <- proto[NA_integer_, , drop = FALSE]
    row$rep <- as.integer(r)
    row$error <- msg
    row
  }
  safe <- function(r) {
    tryCatch({
      out <- one_rep(r)
      if (is.null(out$error)) out$error <- rep(NA_character_, nrow(out))
      out
    }, error = function(e) na_row(r, conditionMessage(e)))
  }
  rows <- if (cores > 1L) {
    parallel::mclapply(reps, safe, mc.cores = cores)
  } else {
    lapply(reps, safe)
  }
  # A crashed worker returns NULL or a try-error instead of a row.
  rows <- Map(function(row, r) {
    if (is.data.frame(row)) row else na_row(r, paste("worker failed:", as.character(row)))
  }, rows, reps)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# The single kinship matrix used by the whole kinship study.
.kinship_G <- function() {
  s <- SETTINGS$kinship
  make_kinship(length(s$rep_counts), n_snp = s$n_snp, n_family = s$n_family,
               seed = BASE_SEED, ridge = s$ridge)
}

# Replicate r of a scenario with one genetic effect (iid, or kinship G).
# h2 = NULL uses the scenario's (first) h2.
.scenario_data <- function(scenario, r, h2 = NULL, G = NULL) {
  s <- SETTINGS[[scenario]]
  stopifnot(!is.null(s$rep_counts))
  if (is.null(h2)) h2 <- s$h2[1]
  seed <- BASE_SEED + r
  design <- make_design(s$rep_counts, seed = seed)
  sim <- simulate_trial(design, h2 = h2, G = G, seed = seed)
  dm <- design_matrices(design)
  list(y = sim$y, g = sim$g, X = dm$X, Z = dm$Z, G = G, h2 = h2)
}

# Predicted genetic values for many gammas at once (one column per gamma).
# Same values as blup(y, X, Z, gamma, G)$g, via
#   g_hat(gamma) = gamma G Z' P_H y,   P_H y = K U (I + gamma D)^-1 U'K'y,
# where K'ZGZ'K = U D U'. One eigen-decomposition serves the whole grid.
.grid_predictions <- function(y, X, Z, G, gammas) {
  K <- qr.Q(qr(X), complete = TRUE)[, -seq_len(ncol(X)), drop = FALSE]
  ZG <- if (is.null(G)) Z else Z %*% G
  eig <- eigen(crossprod(K, tcrossprod(ZG, Z) %*% K), symmetric = TRUE)
  d <- pmax(eig$values, 0)
  KU <- K %*% eig$vectors
  w <- drop(crossprod(KU, y))
  B <- crossprod(ZG, KU)                         # G Z' K U  (m x (n - p))
  (B %*% (w / (1 + outer(d, gammas)))) * rep(gammas, each = nrow(B))
}

# Accuracy and scale of g_hat at every grid point, for one dataset.
# Returns a matrix with rows cor, topk, mse, gain_slope, degenerate.
.curve_metrics <- function(dat, gammas) {
  g_hat <- .grid_predictions(dat$y, dat$X, dat$Z, dat$G, gammas)
  vapply(seq_along(gammas), function(j) {
    gh <- g_hat[, j]
    acc <- accuracy_metrics(gh, dat$g, gh)
    sc <- scale_metrics(gh, dat$g, NA_real_, NA_real_, NA_real_)
    c(cor = acc$cor, topk = acc$topk, mse = mean((gh - dat$g)^2),
      gain_slope = sc$gain_slope, degenerate = acc$degenerate)
  }, numeric(5))
}

# ---------------------------------------------------------------------------
# Experiment: two analysts (hook, Myth 1, cross-validation)
# ---------------------------------------------------------------------------

# 5-fold cross-validation MSE of a model over plots: refit on the training
# folds, predict each held-out plot as X beta + its levels' predicted effects.
# All levels are kept, so a level unseen in training is predicted as 0.
.cv_mse <- function(y, X, effects, folds) {
  pred <- numeric(length(y))
  for (k in sort(unique(folds))) {
    test <- folds == k
    eff_tr <- lapply(effects, function(e) list(Z = e$Z[!test, , drop = FALSE], G = e$G))
    X_tr <- X[!test, , drop = FALSE]
    fit <- reml_dense_fit(y[!test], X_tr, eff_tr)
    b <- blup_general(y[!test], X_tr, eff_tr, fit$gammas)
    p <- drop(X[test, , drop = FALSE] %*% b$beta)
    for (nm in names(effects)) {
      p <- p + drop(effects[[nm]]$Z[test, , drop = FALSE] %*% b$u[[nm]])
    }
    pred[test] <- p
  }
  mean((y - pred)^2)
}

# Offset that gives the CV folds their own random stream. make_design() calls
# set.seed(seed) and starts with sample() over the plots, so reseeding with the
# same seed here would repeat that permutation and tie each plot's fold to its
# genotype's position in the layout, identically in every replicate.
CV_SEED_OFFSET <- 7919L

# Random fold labels (balanced, sizes differ by at most 1) for one replicate.
.cv_folds <- function(seed, n) {
  set.seed(seed + CV_SEED_OFFSET)
  sample(rep_len(seq_len(SETTINGS$cv_folds), n))
}

.two_analysts_rep <- function(truth, r) {
  s <- SETTINGS$analysts
  seed <- BASE_SEED + r
  if (truth == "spurious_family") {
    # No family effect exists; Analyst A fits one anyway.
    design <- make_design(s$rep_counts, family_size = s$family_size, seed = seed)
    sim <- simulate_trial(design, h2 = s$h2, sigma2_family = s$sigma2_family, seed = seed)
    extra <- "family"
    W <- design_matrices(design)$W_family
  } else {
    # A real block effect exists; Analyst A models it.
    design <- make_design(s$rep_counts, block_size = s$block_size, seed = seed)
    sim <- simulate_trial(design, h2 = s$h2, sigma2_block = s$sigma2_block, seed = seed)
    extra <- "block"
    W <- design_matrices(design)$W_block
  }
  dm <- design_matrices(design)
  y <- sim$y
  X <- dm$X

  eff_B <- list(genotype = list(Z = dm$Z, G = NULL))
  eff_A <- c(eff_B, setNames(list(list(Z = W, G = NULL)), extra))
  fit_A <- reml_dense_fit(y, X, eff_A)
  fit_B <- reml_dense_fit(y, X, eff_B)
  u_A <- blup_general(y, X, eff_A, fit_A$gammas)$u
  u_B <- blup_general(y, X, eff_B, fit_B$gammas)$u

  # Predicted genetic value of each genotype (in the order of sim$g).
  pred_B <- unname(u_B$genotype)
  pred_A <- unname(u_A$genotype)
  if (extra == "family") {
    # A believes in families, so a genotype's value includes its family effect.
    fam <- design$family[match(levels(design$genotype), design$genotype)]
    pred_A <- pred_A + unname(u_A$family[as.character(fam)])
  }
  acc_A <- accuracy_metrics(pred_A, sim$g, pred_B)
  acc_B <- accuracy_metrics(pred_B, sim$g, pred_B)

  # Same fold assignment for both analysts.
  folds <- .cv_folds(seed, length(y))

  data.frame(
    rep = as.integer(r),
    loglik_A = fit_A$loglik, loglik_B = fit_B$loglik,
    lrt = 2 * (fit_A$loglik - fit_B$loglik),
    sigma2_extra_hat = unname(fit_A$sigma2[[extra]]),
    cor_A = acc_A$cor, cor_B = acc_B$cor,
    topk_A = acc_A$topk, topk_B = acc_B$topk,
    cv_mse_A = .cv_mse(y, X, eff_A, folds),
    cv_mse_B = .cv_mse(y, X, eff_B, folds),
    degenerate_A = acc_A$degenerate, degenerate_B = acc_B$degenerate
  )
}

# Analyst B fits genotype only; Analyst A adds a family term (truth
# "spurious_family": no family effect exists) or a block term (truth
# "real_block": the block effect is real). Both fits use reml_dense_fit(), so
# their log-likelihoods share the same constants.
run_two_analysts <- function(truth = c("spurious_family", "real_block"), n_rep,
                             cores = 1L) {
  truth <- match.arg(truth)
  proto <- data.frame(
    rep = integer(), loglik_A = numeric(), loglik_B = numeric(), lrt = numeric(),
    sigma2_extra_hat = numeric(), cor_A = numeric(), cor_B = numeric(),
    topk_A = numeric(), topk_B = numeric(), cv_mse_A = numeric(),
    cv_mse_B = numeric(), degenerate_A = logical(), degenerate_B = logical()
  )
  .run_reps(seq_len(n_rep), function(r) .two_analysts_rep(truth, r), proto, cores)
}

# ---------------------------------------------------------------------------
# Experiment: replicate trials of one design (Myth 1b, Myth 3)
# ---------------------------------------------------------------------------

# One replicate: REML fit, plug-in BLUP at gamma_hat, oracle BLUP at the true
# gamma. At the boundary (gamma_hat = 0) log_gamma_hat is -Inf, h2_hat is 0
# and the plug-in predictions are all 0 (degenerate: cor and topk are NA).
.replicate_trial_rep <- function(scenario, r) {
  dat <- .scenario_data(scenario, r)
  setup <- reml_setup(dat$y, dat$X, dat$Z)
  fit <- reml_fit(setup)
  plug <- blup(dat$y, dat$X, dat$Z, fit$gamma)
  oracle <- blup(dat$y, dat$X, dat$Z, dat$h2 / (1 - dat$h2))
  acc <- accuracy_metrics(plug$g, dat$g, oracle$g)
  sc <- scale_metrics(plug$g, dat$g, plug$pev_scaled, fit$sigma2_g, fit$sigma2_e)
  # The curvature is rounding-noise dominated near gamma ~ 1e-4 and can come
  # out <= 0 there; it is only kept where it is a usable (positive) number.
  info <- fit$info_log_gamma
  if (fit$boundary || !isTRUE(info > 0)) info <- NA_real_
  data.frame(
    rep = as.integer(r), loglik_max = fit$loglik, log_gamma_hat = log(fit$gamma),
    h2_hat = fit$h2, boundary = fit$boundary, info_log_gamma = info,
    cor = acc$cor, topk = acc$topk, d_cor = acc$d_cor, d_mse = acc$d_mse,
    gain_slope = sc$gain_slope, reliability = sc$reliability,
    degenerate = acc$degenerate
  )
}

run_replicate_trials <- function(scenario, n_rep, cores = 1L) {
  proto <- data.frame(
    rep = integer(), loglik_max = numeric(), log_gamma_hat = numeric(),
    h2_hat = numeric(), boundary = logical(), info_log_gamma = numeric(),
    cor = numeric(), topk = numeric(), d_cor = numeric(), d_mse = numeric(),
    gain_slope = numeric(), reliability = numeric(), degenerate = logical()
  )
  .run_reps(seq_len(n_rep), function(r) .replicate_trial_rep(scenario, r), proto, cores)
}

# ---------------------------------------------------------------------------
# Experiment: accuracy along the likelihood hill (Myth 2)
# ---------------------------------------------------------------------------

# One dataset (replicate `rep`) evaluated at every gamma of the grid.
# Attributes: gamma_hat (REML estimate, 0 at the boundary), gamma_true.
run_curves <- function(scenario, kinship = FALSE, h2 = NULL, rep = 1) {
  G <- if (kinship) .kinship_G() else NULL
  dat <- .scenario_data(scenario, rep, h2 = h2, G = G)
  setup <- reml_setup(dat$y, dat$X, dat$Z, G)
  fit <- reml_fit(setup)
  grid <- SETTINGS$log_gamma_grid
  gammas <- exp(grid)
  m <- .curve_metrics(dat, gammas)
  out <- data.frame(log_gamma = grid, h2 = gammas / (1 + gammas),
                    loglik = reml_loglik(gammas, setup),
                    cor = m["cor", ], topk = m["topk", ], mse = m["mse", ],
                    gain_slope = m["gain_slope", ])
  attr(out, "gamma_hat") <- fit$gamma
  attr(out, "gamma_true") <- dat$h2 / (1 - dat$h2)
  out
}

# The same grid averaged over replicates 1..n_rep. cor, topk and gain_slope
# average the non-degenerate replicates at each grid point (n_valid of them;
# NA when none); mse averages every successful replicate. Failed replicates
# are excluded and reported in the attributes n_failed and errors.
run_expected_curves <- function(scenario, kinship = FALSE, h2 = NULL, n_rep,
                                cores = 1L) {
  G <- if (kinship) .kinship_G() else NULL
  grid <- SETTINGS$log_gamma_grid
  gammas <- exp(grid)
  one <- function(r) {
    tryCatch(.curve_metrics(.scenario_data(scenario, r, h2 = h2, G = G), gammas),
             error = function(e) conditionMessage(e))
  }
  per_rep <- if (cores > 1L) {
    parallel::mclapply(seq_len(n_rep), one, mc.cores = cores)
  } else {
    lapply(seq_len(n_rep), one)
  }
  ok <- vapply(per_rep, is.matrix, logical(1))
  errors <- vapply(per_rep[!ok], function(x) paste(as.character(x), collapse = " "),
                   character(1))
  per_rep <- per_rep[ok]

  stack <- function(name) vapply(per_rep, function(m) m[name, ], numeric(length(gammas)))
  as_matrix <- function(x) matrix(x, nrow = length(gammas))
  valid <- as_matrix(stack("degenerate")) == 0
  n_valid <- rowSums(valid)
  mean_valid <- function(name) {
    x <- as_matrix(stack(name))
    x[!valid] <- NA
    ifelse(n_valid > 0, rowSums(x, na.rm = TRUE) / n_valid, NA_real_)
  }
  out <- data.frame(log_gamma = grid, h2 = gammas / (1 + gammas),
                    cor = mean_valid("cor"), topk = mean_valid("topk"),
                    mse = if (any(ok)) rowMeans(as_matrix(stack("mse"))) else NA_real_,
                    gain_slope = mean_valid("gain_slope"),
                    n_valid = as.integer(n_valid))
  h2_true <- if (is.null(h2)) SETTINGS[[scenario]]$h2[1] else h2
  attr(out, "gamma_true") <- h2_true / (1 - h2_true)
  attr(out, "n_rep") <- as.integer(n_rep)
  attr(out, "n_failed") <- sum(!ok)
  attr(out, "errors") <- errors
  out
}

# ---------------------------------------------------------------------------
# Experiment: convergence (Myth 4)
# ---------------------------------------------------------------------------

# AI-REML on the baseline scenario (one genetic component) from each start and
# iteration cap. loglik is the solver's own (unprofiled) log-likelihood at the
# point it returns; loglik_gap = loglik - REML maximum (<= 0 up to rounding).
run_convergence <- function(n_rep, maxit = c(1, 2, 3, 30),
                            starts = list(default = NULL,
                                          poor = c(genotype = 0.001, residual = 10)),
                            cores = 1L) {
  proto <- data.frame(
    rep = integer(), start = character(), maxit = integer(),
    converged = logical(), loglik = numeric(), sigma2_g = numeric(),
    sigma2_e = numeric(), h2 = numeric(), cor = numeric(), loglik_gap = numeric(),
    error = character()
  )
  cells <- expand.grid(maxit = as.integer(maxit), start = names(starts),
                       stringsAsFactors = FALSE)
  one_rep <- function(r) {
    dat <- .scenario_data("baseline", r)
    best <- reml_fit(reml_setup(dat$y, dat$X, dat$Z))$loglik
    effects <- list(genotype = list(Z = dat$Z, G = NULL))
    rows <- lapply(seq_len(nrow(cells)), function(i) {
      st <- cells$start[i]
      mi <- cells$maxit[i]
      tryCatch({
        # Only the two expected non-convergence warnings (small maxit) are muffled.
        fit <- withCallingHandlers(
          ai_reml(dat$y, dat$X, effects, start = starts[[st]], maxit = mi),
          warning = function(w) {
            if (grepl("not converged|changed by more than", conditionMessage(w))) {
              invokeRestart("muffleWarning")
            }
          })
        s2g <- fit$sigma2[["genotype"]]
        s2e <- fit$sigma2[["residual"]]
        g_hat <- blup(dat$y, dat$X, dat$Z, s2g / s2e)$g
        data.frame(rep = as.integer(r), start = st, maxit = mi,
                   converged = fit$converged, loglik = fit$loglik,
                   sigma2_g = s2g, sigma2_e = s2e, h2 = s2g / (s2g + s2e),
                   cor = accuracy_metrics(g_hat, dat$g, g_hat)$cor,
                   loglik_gap = fit$loglik - best, error = NA_character_)
      }, error = function(e) {
        row <- proto[NA_integer_, , drop = FALSE]
        row$rep <- as.integer(r)
        row$start <- st
        row$maxit <- mi
        row$error <- conditionMessage(e)
        row
      })
    })
    do.call(rbind, rows)
  }
  .run_reps(seq_len(n_rep), one_rep, proto[setdiff(names(proto), "error")], cores)
}
