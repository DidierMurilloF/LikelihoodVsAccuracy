# Dense REML oracle and BLUP for one random genetic effect.
#
# Model:  y = X beta + Z g + e,   g ~ N(0, sigma2_g G),   e ~ N(0, sigma2_e I)
# gamma = sigma2_g / sigma2_e (signal-to-noise), lambda = 1 / gamma.
# G = NULL means G = I (independent genotypes).

# One-time preparation for a dataset: everything that does not depend on gamma.
reml_setup <- function(y, X, Z, G = NULL) {
  n <- length(y)
  p <- ncol(X)
  qx <- qr(X)
  stopifnot(qx$rank == p)
  # K: orthonormal contrasts free of the fixed effects (K'X = 0, K'K = I).
  K <- qr.Q(qx, complete = TRUE)[, -seq_len(p), drop = FALSE]
  ZGZ <- if (is.null(G)) tcrossprod(Z) else Z %*% G %*% t(Z)
  # Eigen-decomposition K'ZGZ'K = U D U'.
  eig <- eigen(crossprod(K, ZGZ %*% K), symmetric = TRUE)
  d <- pmax(eig$values, 0)
  w <- drop(crossprod(eig$vectors, crossprod(K, y)))
  list(d = d, w = w, n = n, p = p, y = y, X = X, Z = Z, G = G)
}

# REML estimate of sigma2_e for a given gamma (vectorized over gamma).
reml_sigma2_e <- function(gamma, setup) {
  vapply(gamma, function(g) {
    sum(setup$w^2 / (1 + g * setup$d)) / (setup$n - setup$p)
  }, numeric(1))
}

# Profiled REML log-likelihood as a function of gamma (vectorized over gamma).
reml_loglik <- function(gamma, setup) {
  nu <- setup$n - setup$p                        # residual degrees of freedom
  s2e <- reml_sigma2_e(gamma, setup)             # best sigma2_e for each gamma
  logdet <- vapply(gamma, function(g) sum(log1p(g * setup$d)), numeric(1))
  -0.5 * (nu * (log(2 * pi) + 1) + nu * log(s2e) + logdet)
}

# Find the REML estimate: grid search on log(gamma) plus gamma = 0,
# then optimize() between the neighbours of the best grid point.
reml_fit <- function(setup, log_gamma_grid = SETTINGS$log_gamma_grid) {
  ll_fun <- function(lg) reml_loglik(exp(lg), setup)
  grid_ll <- ll_fun(log_gamma_grid)
  best <- which.max(grid_ll)
  bracket <- log_gamma_grid[c(max(best - 1, 1), min(best + 1, length(log_gamma_grid)))]
  opt <- optimize(ll_fun, bracket, maximum = TRUE, tol = 1e-10)

  # Keep the better of the grid point and the optimizer.
  if (opt$objective >= grid_ll[best]) {
    log_gamma <- opt$maximum
    loglik <- opt$objective
  } else {
    log_gamma <- log_gamma_grid[best]
    loglik <- grid_ll[best]
  }
  gamma <- exp(log_gamma)

  # Boundary: no genetic variance (gamma = 0) beats every interior point.
  ll_zero <- reml_loglik(0, setup)
  boundary <- ll_zero > loglik
  if (boundary) {
    gamma <- 0
    loglik <- ll_zero
  }

  # Curvature in log(gamma): how sharply the data pin gamma down.
  h <- 1e-3
  info <- if (boundary) NA_real_ else
    -(ll_fun(log_gamma + h) - 2 * loglik + ll_fun(log_gamma - h)) / h^2

  sigma2_e <- reml_sigma2_e(gamma, setup)
  list(gamma = gamma, sigma2_g = gamma * sigma2_e, sigma2_e = sigma2_e,
       h2 = gamma / (1 + gamma), loglik = loglik, boundary = boundary,
       info_log_gamma = info,
       grid = data.frame(log_gamma = log_gamma_grid, loglik = grid_ll))
}

# Reference implementation: the same REML log-likelihood from the full n x n
# matrix V (Harville form R1, plus 0.5 log|X'X| to match the contrast form).
reml_loglik_naive <- function(sigma2_g, sigma2_e, y, X, Z, G = NULL) {
  n <- length(y)
  p <- ncol(X)
  ZGZ <- if (is.null(G)) tcrossprod(Z) else Z %*% G %*% t(Z)
  V <- sigma2_g * ZGZ + sigma2_e * diag(n)
  R <- chol(V)
  Vi <- chol2inv(R)
  XViX <- crossprod(X, Vi %*% X)
  P <- Vi - Vi %*% X %*% solve(XViX, crossprod(X, Vi))
  logdet <- function(M) as.numeric(determinant(M, logarithm = TRUE)$modulus)
  -0.5 * ((n - p) * log(2 * pi) + 2 * sum(log(diag(R))) + logdet(XViX) +
            drop(crossprod(y, P %*% y))) + 0.5 * logdet(crossprod(X))
}

# BLUP of genetic values from Henderson's mixed model equations (MME).
blup <- function(y, X, Z, gamma, G = NULL) {
  p <- ncol(X)
  m <- ncol(Z)
  # gamma = 0: no genetic variance, so every BLUP is 0 and beta is least squares.
  if (gamma == 0) {
    beta <- drop(solve(crossprod(X), crossprod(X, y)))
    return(list(beta = beta, g = rep(0, m), pev_scaled = rep(NA_real_, m)))
  }
  lambda <- 1 / gamma
  Ginv <- if (is.null(G)) diag(m) else solve(G)
  # MME:  [X'X  X'Z             ] [beta]   [X'y]
  #       [Z'X  Z'Z + lambda G^-1] [g   ] = [Z'y]
  C <- rbind(cbind(crossprod(X),    crossprod(X, Z)),
             cbind(crossprod(Z, X), crossprod(Z) + lambda * Ginv))
  rhs <- c(crossprod(X, y), crossprod(Z, y))
  C_inv <- solve(C)
  sol <- drop(C_inv %*% rhs)
  g_idx <- p + seq_len(m)
  # Diagonal of the g-block of C^-1 = prediction error variance / sigma2_e.
  list(beta = sol[seq_len(p)], g = sol[g_idx], pev_scaled = diag(C_inv)[g_idx])
}
