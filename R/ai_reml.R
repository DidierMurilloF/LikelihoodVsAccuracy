# Average-information REML (AI-REML) in the variance components themselves,
# with ASReml-style convergence warnings.
#
# Model:  y = X beta + sum_k Z_k u_k + e,
#         u_k ~ N(0, sigma2_k G_k),   e ~ N(0, sigma2_e I)
# theta = (sigma2_1, ..., sigma2_K, sigma2_e) and
#         V = sum_k sigma2_k V_k + sigma2_e I,   V_k = Z_k G_k Z_k'.
# `effects` follows R/oracle_multi.R: a named list of list(Z = , G = ), with
# G = NULL meaning G = I.
#
# Everything is computed with error contrasts K'y (K'X = 0, K'K = I), so that
# P = K (K'VK)^-1 K' and the log-likelihood is the contrast form
#   l_K(theta) = -1/2 [ (n-p) log 2pi + log|K'VK| + y'K (K'VK)^-1 K'y ],
# with the same constants as reml_loglik() and reml_loglik_naive().
# Unlike reml_loglik_multi(), sigma2_e is NOT profiled out: l_K is evaluated
# at the point theta exactly as given, so an unconverged fit reports its own,
# honest log-likelihood.

# One-time preparation: contrasts K'y and A_k = K' V_k K for every component,
# with the residual (V_e = I, so A_e = I) as the last one.
.ai_setup <- function(y, X, effects) {
  stopifnot(is.list(effects), length(effects) >= 1,
            !is.null(names(effects)), all(nzchar(names(effects))),
            !anyDuplicated(names(effects)), !"residual" %in% names(effects))
  p <- ncol(X)
  qx <- qr(X)
  stopifnot(qx$rank == p)
  K <- qr.Q(qx, complete = TRUE)[, -seq_len(p), drop = FALSE]
  nu <- length(y) - p
  A <- lapply(effects, function(eff) crossprod(K, .zgz(eff) %*% K))
  A$residual <- diag(nu)
  list(w = drop(crossprod(K, y)), A = A, nu = nu)
}

# Check a named theta against the components and return it in their order.
.ai_theta <- function(sigma2, setup) {
  nms <- names(setup$A)
  stopifnot(is.numeric(sigma2), setequal(names(sigma2), nms))
  sigma2[nms]
}

# Log-likelihood, score and average-information matrix at theta.
# With C = K'VK and a = C^-1 K'y (so that P y = K a):
#   score  s_k   = -1/2 [ tr(P V_k) - y'P V_k P y ]
#                = -1/2 [ tr(C^-1 A_k) - a'A_k a ]
#   AI     AI_kl =  1/2 y'P V_k P V_l P y = 1/2 (A_k a)' C^-1 (A_l a)
.ai_eval <- function(theta, setup) {
  C <- Reduce(`+`, Map(`*`, theta, setup$A))
  R <- chol(C)
  Cinv <- function(M) backsolve(R, backsolve(R, M, transpose = TRUE))
  a <- Cinv(setup$w)
  loglik <- -0.5 * (setup$nu * log(2 * pi) + 2 * sum(log(diag(R))) +
                      sum(setup$w * a))

  Ci <- chol2inv(R)
  Aa <- vapply(setup$A, function(Ak) drop(Ak %*% a), numeric(setup$nu))
  score <- -0.5 * (vapply(setup$A, function(Ak) sum(Ci * Ak), numeric(1)) -
                     colSums(a * Aa))
  AI <- 0.5 * crossprod(Aa, Cinv(Aa))
  list(loglik = loglik, score = score, AI = AI)
}

# Unprofiled REML log-likelihood l_K at sigma2 (named, including `residual`).
reml_loglik_theta <- function(sigma2, y, X, effects) {
  setup <- .ai_setup(y, X, effects)
  .ai_eval(.ai_theta(sigma2, setup), setup)$loglik
}

# REML score (gradient of l_K) with respect to sigma2 (named, incl. `residual`).
reml_score_theta <- function(sigma2, y, X, effects) {
  setup <- .ai_setup(y, X, effects)
  .ai_eval(.ai_theta(sigma2, setup), setup)$score
}

# AI-REML fit.
#   start     named starting values (same names as the returned sigma2);
#             default var(y) / (k + 1) for each of the k components and the
#             residual.
#   maxit     maximum number of AI updates.
#   tol_loglik, tol_par
#             converged when |change in loglik| < tol_loglik AND every
#             relative change |delta theta_k| / theta_k < tol_par. Components
#             held at the lower bound are left out of the relative change
#             (they have converged at the boundary).
# Returns list(sigma2, loglik, iterations, converged, history), where history
# has one row per iteration (iter 0 = the start) and loglik is l_K at the
# returned sigma2. If maxit is reached without convergence it warns like
# ASReml does.
ai_reml <- function(y, X, effects, start = NULL, maxit = 30,
                    tol_loglik = 1e-6, tol_par = 0.01) {
  setup <- .ai_setup(y, X, effects)
  nms <- names(setup$A)
  floor_val <- 1e-8 * var(y)           # lower bound for every component
  theta <- if (is.null(start)) {
    setNames(rep(var(y) / length(nms), length(nms)), nms)
  } else {
    .ai_theta(start, setup)
  }
  stopifnot(all(theta > 0))
  theta <- pmax(theta, floor_val)

  ev <- .ai_eval(theta, setup)
  history <- list(c(iter = 0, loglik = ev$loglik, theta))
  converged <- FALSE
  par_ok <- FALSE
  iter <- 0

  while (iter < maxit && !converged) {
    iter <- iter + 1

    # Newton-type step theta + AI^-1 s. A component already at the floor that
    # the step would push further down is held there, and the step is
    # recomputed for the remaining (free) components.
    free <- rep(TRUE, length(theta))
    repeat {
      delta <- rep(0, length(theta))
      delta[free] <- solve(ev$AI[free, free, drop = FALSE], ev$score[free])
      stuck <- free & theta <= floor_val & delta < 0
      if (!any(stuck)) break
      free <- free & !stuck
    }

    # Step-halving: shorten the step while any component would fall below the
    # floor; after 10 halvings, clamp at the floor.
    step <- 1
    for (h in seq_len(10)) {
      if (all(theta + step * delta >= floor_val)) break
      step <- step / 2
    }
    theta_new <- pmax(theta + step * delta, floor_val)
    ev_new <- .ai_eval(theta_new, setup)

    # Convergence: log-likelihood change and relative parameter changes.
    interior <- theta_new > floor_val
    rel_change <- abs(theta_new - theta)[interior] / theta_new[interior]
    loglik_ok <- abs(ev_new$loglik - ev$loglik) < tol_loglik
    par_ok <- all(rel_change < tol_par)
    converged <- loglik_ok && par_ok

    theta <- theta_new
    ev <- ev_new
    history[[iter + 1]] <- c(iter = iter, loglik = ev$loglik, theta)
  }

  if (!converged) {
    warning("Log-likelihood not converged", call. = FALSE)
    if (!par_ok) {
      warning(sprintf("Some components changed by more than %s%% on the last iteration",
                      format(100 * tol_par)), call. = FALSE)
    }
  }

  history <- as.data.frame(do.call(rbind, history))
  history$iter <- as.integer(history$iter)
  list(sigma2 = theta, loglik = ev$loglik, iterations = iter,
       converged = converged, history = history)
}
