# Simulation helpers: field designs, VanRaden kinship, and trial data.
# Scale: sigma2_g = h2, sigma2_e = 1 - h2, mu = 0 (estimated later as X = 1).

# Field design. Rows are plots in random field order (plot = 1..n).
# Families group genotypes (genotype i is in family ceiling(i / family_size));
# blocks are consecutive runs of `block_size` plots in field order.
make_design <- function(rep_counts, family_size = NULL, block_size = NULL, seed) {
  set.seed(seed)
  m <- length(rep_counts)
  n <- sum(rep_counts)
  geno_levels <- sprintf("G%03d", seq_len(m))
  geno_id <- sample(rep(seq_len(m), rep_counts))
  genotype <- factor(geno_levels[geno_id], levels = geno_levels)

  if (is.null(family_size)) {
    family <- factor(rep(NA_character_, n))
  } else {
    n_family <- ceiling(m / family_size)
    fam_levels <- sprintf("F%02d", seq_len(n_family))
    family <- factor(fam_levels[ceiling(geno_id / family_size)], levels = fam_levels)
  }

  if (is.null(block_size)) {
    block <- factor(rep(NA_character_, n))
  } else {
    n_block <- ceiling(n / block_size)
    blk_levels <- sprintf("B%03d", seq_len(n_block))
    block <- factor(blk_levels[ceiling(seq_len(n) / block_size)], levels = blk_levels)
  }

  data.frame(plot = seq_len(n), genotype = genotype, family = family, block = block)
}

# VanRaden (method 1) kinship from simulated SNPs with family-structured allele
# frequencies, rescaled to mean diagonal 1 (so Var(g_i) = sigma2_g on average),
# plus ridge * I added after rescaling. Genotypes are split evenly into
# n_family families.
make_kinship <- function(m, n_snp = 1000, n_family = 10, seed, ridge = 1e-3) {
  set.seed(seed)
  fam <- ceiling(seq_len(m) * n_family / m)
  M <- matrix(0, m, n_snp)
  for (f in seq_len(n_family)) {
    idx <- which(fam == f)
    p_f <- runif(n_snp, 0.05, 0.95)
    M[idx, ] <- matrix(rbinom(length(idx) * n_snp, 2, rep(p_f, each = length(idx))),
                       length(idx), n_snp)
  }
  p <- colMeans(M) / 2
  keep <- p > 0 & p < 1               # drop monomorphic markers
  M <- M[, keep, drop = FALSE]
  p <- p[keep]
  W <- sweep(M, 2, 2 * p)
  G <- tcrossprod(W) / (2 * sum(p * (1 - p)))
  G <- (G + t(G)) / 2
  G <- G / mean(diag(G))
  G + ridge * diag(m)
}

# Indicator matrix of a factor with all levels kept; NULL if absent / all NA.
.incidence <- function(f) {
  if (is.null(f) || all(is.na(f))) return(NULL)
  f <- factor(f)
  out <- matrix(0, length(f), nlevels(f))
  out[cbind(seq_along(f), as.integer(f))] <- 1
  colnames(out) <- levels(f)
  out
}

design_matrices <- function(design) {
  geno <- factor(design$genotype)
  Z <- .incidence(geno)
  list(
    X = matrix(1, nrow(design), 1),
    Z = Z,
    W_family = .incidence(design$family),
    W_block = .incidence(design$block)
  )
}

# Simulate one trial. Family / block effects are drawn only when their variance
# is positive, so a zero variance leaves y identical to a design without them.
simulate_trial <- function(design, h2, G = NULL, sigma2_family = 0,
                           sigma2_block = 0, seed) {
  set.seed(seed)
  geno <- factor(design$genotype)
  m <- nlevels(geno)
  n <- nrow(design)
  s2g <- h2
  s2e <- 1 - h2

  if (is.null(G)) {
    g <- rnorm(m, 0, sqrt(s2g))
  } else {
    L <- t(chol(G))
    g <- sqrt(s2g) * drop(L %*% rnorm(m))
  }
  y <- g[as.integer(geno)]

  if (sigma2_family > 0) {
    fam <- factor(design$family)
    u <- rnorm(nlevels(fam), 0, sqrt(sigma2_family))
    y <- y + u[as.integer(fam)]
  }
  if (sigma2_block > 0) {
    blk <- factor(design$block)
    v <- rnorm(nlevels(blk), 0, sqrt(sigma2_block))
    y <- y + v[as.integer(blk)]
  }
  y <- y + rnorm(n, 0, sqrt(s2e))

  list(y = y, g = g, design = design)
}
