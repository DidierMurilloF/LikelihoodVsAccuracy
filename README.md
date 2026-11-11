# Likelihood vs Prediction Accuracy in Mixed Models

A higher REML log-likelihood does not guarantee more accurate predictions of genotype effects. This article demonstrates that with simulation, using a transparent base-R oracle (exact BLUP and truth-based accuracy) and a from-scratch AI-REML implementation, so every number can be traced and reproduced without ASReml, lme4, or similar packages.

## Reproduce

Built with R 4.6.2. The base-R code in `R/` needs no packages; the article
and the banner need:

```r
install.packages(c("distill", "rmarkdown", "dplyr", "tidyr", "ggplot2",
                   "patchwork", "kableExtra", "plotly"))
```

1. Run the simulations (writes to `results/`; about 1 hour on 2 cores):

   ```sh
   Rscript scripts/run_simulations.R
   ```

   The replicates run in parallel on `LVA_CORES` cores (default 2, capped at
   the number of cores minus one), for example
   `LVA_CORES=4 Rscript scripts/run_simulations.R`. Every replicate sets its
   own seed, so the results are identical for any core count.

2. Compute the profile-likelihood intervals used in Myth 3 (writes
   `results/profile_coverage.rds`; same seeds as step 1, about 2.5 minutes on
   2 cores; also honours `LVA_CORES`):

   ```sh
   Rscript scripts/profile_coverage.R
   ```

3. Run the robustness checks used in Myth 2 (uneven replication, and
   within- vs between-family accuracy with the genomic relationship matrix;
   writes `results/robustness_checks.rds`; 150 trials per setting, about
   1 minute on 2 cores; also honours `LVA_CORES`, and `N_REP` sets the
   number of trials):

   ```sh
   Rscript scripts/robustness_checks.R
   ```

## License

The code (`R/`, `scripts/`) is released under the MIT License (see
`LICENSE`). The article text and figures are licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
