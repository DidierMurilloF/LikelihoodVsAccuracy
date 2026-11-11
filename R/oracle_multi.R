# Dense REML and BLUP for several independent random effects.
#
# Model:  y = X beta + sum_k Z_k u_k + e,
#         u_k ~ N(0, sigma2_k G_k),   e ~ N(0, sigma2_e I)
# gamma_k = sigma2_k / sigma2_e >= 0, and V = sigma2_e H with
#         H = I + sum_k gamma_k Z_k G_k Z_k'.
#
# `effects` is a named list; each element is list(Z = <n x q>, G = <q x q or
# NULL>), with G = NULL meaning G = I. The genotype effect is named "genotype".

# Z G Z' for one effect.
.zgz <- function(eff) {
  if (is.null(eff$G)) tcrossprod(eff$Z) else eff$Z %*% eff$G %*% t(eff$Z)
}

# One-time preparation for a dataset: everything that does not depend on the
# gammas. Works with error contrasts K'y (K'X = 0, K'K = I), so that
#   log|K'HK| = log|H| + log|X'H^-1 X| - log|X'X|   and   P_H = K (K'HK)^-1 K'.
.multi_setup <- function(y, X, effects) {
  stopifnot(is.list(effects), length(effects) >= 1,
            !is.null(names(effects)), all(nzchar(names(effects))),
            !anyDuplicated(names(effects)))
  n <- length(y)
  p <- ncol(X)
  qx <- qr(X)
  stopifnot(qx$rank == p)
  K <- qr.Q(qx, complete = TRUE)[, -seq_len(p), drop = FALSE]
  # A_k = K' Z_k G_k Z_k' K, one (n-p) x (n-p) matrix per effect.
  A <- lapply(effects, function(eff) crossprod(K, .zgz(eff) %*% K))
  list(w = drop(crossprod(K, y)), A = A, n = n, p = p)
}

# Profiled REML log-likelihood (contrast form l_K, same constants as
# reml_loglik) and, optionally, its gradient with respect to the gammas:
#   d l / d gamma_k = -1/2 [ tr(C^-1 A_k) - (n-p) a'A_k a / (w'a) ],
# where C = K'HK = I + sum_k gamma_k A_k and a = C^-1 w.
.multi_eval <- function(gammas, setup, gradient = FALSE) {
  nu <- setup$n - setup$p
  C <- diag(nu)
  for (k in seq_along(setup$A)) {
    if (gammas[k] != 0) C <- C + gammas[k] * setup$A[[k]]
  }
  R <- chol(C)
  a <- backsolve(R, backsolve(R, setup$w, transpose = TRUE))
  ywy <- sum(setup$w * a)                       # y'P_H y
  s2e <- ywy / nu
  logdet <- 2 * sum(log(diag(R)))               # log|K'HK|
  loglik <- -0.5 * (nu * (log(2 * pi) + 1) + nu * log(s2e) + logdet)
  if (!gradient) return(list(loglik = loglik, sigma2_e = s2e))

  Ci <- chol2inv(R)
  grad <- vapply(setup$A, function(Ak) {
    -0.5 * (sum(Ci * Ak) - nu * sum(a * (Ak %*% a)) / ywy)
  }, numeric(1))
  list(loglik = loglik, sigma2_e = s2e, gradient = grad)
}

# Profiled REML log-likelihood at gammas = sigma2_k / sigma2_e (in the order
# of `effects`). Same quantity and constants as reml_loglik().
reml_loglik_multi <- function(gammas, y, X, effects) {
  stopifnot(length(gammas) == length(effects), all(gammas >= 0))
  .multi_eval(gammas, .multi_setup(y, X, effects))$loglik
}

# REML fit: L-BFGS-B on the gammas (lower bound 0) from every combination of
# the starts {0, 0.1, 1, 10} per effect; the best result over starts is kept.
reml_dense_fit <- function(y, X, effects) {
  setup <- .multi_setup(y, X, effects)
  nms <- names(effects)
  n_eff <- length(effects)

  # Cache the last evaluation: optim() calls fn and gr at the same point.
  cache <- new.env()
  eval_at <- function(gm) {
    if (!identical(gm, cache$gm)) {
      cache$gm <- gm
      cache$val <- .multi_eval(gm, setup, gradient = TRUE)
    }
    cache$val
  }
  fn <- function(gm) -eval_at(gm)$loglik
  gr <- function(gm) -eval_at(gm)$gradient

  starts <- as.matrix(expand.grid(rep(list(c(0, 0.1, 1, 10)), n_eff)))
  best <- NULL
  for (i in seq_len(nrow(starts))) {
    opt <- optim(unname(starts[i, ]), fn, gr, method = "L-BFGS-B",
                 lower = 0, control = list(factr = 10, maxit = 500))
    if (is.null(best) || -opt$value > best$loglik) {
      best <- list(gammas = pmax(opt$par, 0), loglik = -opt$value)
    }
  }

  # Snap tiny gammas to exactly 0 when that does not lower the log-likelihood.
  gammas <- best$gammas
  loglik <- .multi_eval(gammas, setup)$loglik
  for (k in which(gammas > 0 & gammas < 1e-8)) {
    trial <- gammas
    trial[k] <- 0
    ll_trial <- .multi_eval(trial, setup)$loglik
    if (ll_trial >= loglik) {
      gammas <- trial
      loglik <- ll_trial
    }
  }

  s2e <- .multi_eval(gammas, setup)$sigma2_e
  names(gammas) <- nms
  list(gammas = gammas,
       sigma2 = c(gammas * s2e, residual = s2e),
       loglik = loglik,
       boundary = setNames(gammas == 0, nms))
}

# BLUP of every random effect for known gammas:
#   beta = (X'H^-1 X)^-1 X'H^-1 y (GLS),  u_k = gamma_k G_k Z_k' H^-1 (y - X beta).
blup_general <- function(y, X, effects, gammas) {
  stopifnot(length(gammas) == length(effects), all(gammas >= 0))
  H <- diag(length(y))
  for (k in seq_along(effects)) {
    if (gammas[k] != 0) H <- H + gammas[k] * .zgz(effects[[k]])
  }
  R <- chol(H)
  Hinv <- function(M) backsolve(R, backsolve(R, M, transpose = TRUE))
  HiX <- Hinv(X)
  beta <- drop(solve(crossprod(X, HiX), crossprod(HiX, y)))
  r <- drop(Hinv(y - X %*% beta))               # H^-1 (y - X beta)
  u <- lapply(seq_along(effects), function(k) {
    eff <- effects[[k]]
    Zr <- drop(crossprod(eff$Z, r))
    gammas[k] * (if (is.null(eff$G)) Zr else drop(eff$G %*% Zr))
  })
  names(u) <- names(effects)
  list(beta = beta, u = u)
}
