# Central settings for every experiment in the article.
# All scenario designs are ordered by genotype index (genotype 1 first).

BASE_SEED <- 20261109

.rep_baseline <- c(rep(2L, 50), rep(1L, 150))   # n = 250 plots

SETTINGS <- list(
  # Weak signal: same design as baseline, low heritability.
  weak = list(
    rep_counts = .rep_baseline,
    h2 = 0.1
  ),
  # Baseline: 50 genotypes with 2 plots, 150 with 1 plot.
  baseline = list(
    rep_counts = .rep_baseline,
    h2 = 0.3
  ),
  # Strong signal and a larger design: 100 x 3 plots, 100 x 2 plots (n = 500).
  strong = list(
    rep_counts = c(rep(3L, 100), rep(2L, 100)),
    h2 = 0.5
  ),
  # Analysts' models: spurious family term (truth 0) vs real block term.
  analysts = list(
    rep_counts = c(rep(2L, 40), rep(1L, 160)),   # n = 240 plots
    h2 = 0.3,
    family_size = 10L,       # 20 families x 10 genotypes
    block_size = 6L,         # incomplete blocks of 6 consecutive plots
    sigma2_family = 0,       # spurious-family truth
    sigma2_block = 0.15      # real-block truth
  ),
  # Kinship: baseline design, VanRaden G from SNPs, ridge for invertibility.
  kinship = list(
    rep_counts = .rep_baseline,
    h2 = c(0.3, 0.1),        # first value is the default
    n_snp = 1000L,
    n_family = 10L,          # 10 families x 20 genotypes
    ridge = 1e-3
  ),
  # Shared across scenarios.
  n_rep = 500L,
  top_frac = 0.1,            # top-k with k = 10% of genotypes (20)
  cv_folds = 5L,
  log_gamma_grid = seq(log(1e-4), log(1e4), length.out = 401)
)

rm(.rep_baseline)
