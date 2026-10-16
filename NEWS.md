# iDESC 0.4.0

* Recover coefficient inference for `zi_model = "constant"` fits without a positive definite
  Hessian or group standard errors. At the pi = 0 boundary, use the nested NB fit (when it has
  usable standard errors) for
  coefficients, Wald tests, pairwise contrasts, dispersion and subject variance; retain the
  original ZINB `Theta`. New `zi_boundary_tol` defaults to 1e-3 (twice the absolute likelihood
  difference). Otherwise retry ZINB from NB starting values, then with BFGS, accepting only
  usable fits with a likelihood sufficiently close to the original.
* Add `zi_status` to `gene_status` and immediately after `Status` in `model_results`:
  `none`, `estimated`, `boundary`, `refit`, `unidentified`, or `NA` for genes not fitted.
  `Status` reflects the chosen coefficient inference; convergence/Hessian diagnostics still
  describe the original fit.
* LRTs, likelihood summaries, Npar, AIC and all permutation statistics retain the original
  ZINB fit. Permuted fits never run recovery. NB (`zi_model = "none"`) results are unchanged
  apart from the new status column.

# iDESC 0.3.0

**Breaking release.** The model fitted by `iDESC()` changed, so results differ from 0.2.0
*by design*, and several arguments and output columns were removed or renamed. Version 0.2.0
(commit `2be58b9`) remains available to reproduce earlier analyses:
`remotes::install_github("pavillos/iDESC", ref = "v0.2.0")`.

## Removed

* **The original zero-inflation model**
  `ziformula = ~ offset(predict_pi_offset) - 1 + ind_zero + (1 | sub)` is gone. Its covariate
  `ind_zero = 1 * (gene == 0)` is the outcome used as its own predictor (a cell "predicts" its
  structural zero from the fact that its count is zero), so the model is ill-defined: the
  coefficient of `ind_zero` diverges, the Hessian is singular and most genes end without
  standard errors. On real data it gave Wald p-values for about 20% of the genes only.
* Everything that existed only for that model:
  * the argument `use_ind_zero` (added after 0.2.0, never released);
  * the LOESS prediction of the zero proportion inside `iDESC()` and its arguments `span` and
    `loess_control` (`iDESC()` no longer draws a plot);
  * the exported function `zp_prediction()`. It only computed the offset of the removed model,
    nothing else used it, and it was the only reason for the dependencies on `Seurat` and
    `graphics`; it was removed rather than kept as an orphan.
* Output columns `Sigma2_ZI` (variance of the zero-inflation random effect) and `Deviance`.
* Dependencies `Seurat`, `graphics`, `plyr` and `boot` (no longer used).

## New

* `zi_model = c("none", "constant")`. `"none"` (default) is the negative binomial mixed model
  `y ~ group + covariates + offset(log(norm_sf)) + (1 | sub)` without zero inflation
  (`ziformula = ~ 0`); `"constant"` adds one zero-inflation intercept per gene (`ziformula = ~ 1`,
  reported in `Theta`).
* **Likelihood-ratio test of the group effect**, always computed: the full model against the same
  model without `group`. New columns `LRT_stat`, `Pval_LRT_group` and `FDR_LRT_group` (BH across
  the genes with a p-value). The log-likelihoods are taken from the optimizer
  (`-fit$fit$objective`), so the test is available when the Hessian is not positive definite.
  The Wald test of every coefficient is kept (`Beta_*`, `Pval_Beta_*`, `FDR_Beta_*`).
* **Pairwise contrasts**: `pairwise = TRUE` (default; used when the group has three or more
  levels) adds the contrasts between non-reference levels, `Beta_<B>_vs_<A>`,
  `Pval_Beta_<B>_vs_<A>` (Wald, from the covariance of the fixed effects) and
  `FDR_Beta_<B>_vs_<A>`.
* **Permutation p-values** for designs with few subjects, where the asymptotic p-values are
  anti-conservative. Arguments `permutation = c("auto", "yes", "no")`,
  `min_subjects_asymptotic = 10`, `n_permutations = 999`, `seed = 1`, `sequential = TRUE`.
  * `"auto"` uses the permutation test when the smallest group has fewer than
    `min_subjects_asymptotic` subjects.
  * The group label is permuted between subjects (each subject keeps its cells and covariates,
    group sizes are kept). The null model is fitted once per gene; each permutation refits the
    full model only. The same permutations, generated from `seed`, are used for every gene:
    results are reproducible and independent of `cores`.
  * New columns `Pval_perm_group` / `FDR_perm_group` (likelihood-ratio statistic),
    `Pval_perm_Beta_<contrast>` / `FDR_perm_Beta_<contrast>` (absolute estimate of every group
    coefficient and pairwise contrast) and `N_perm` (valid permutations used for the gene).
    p = (1 + number of permutations with statistic >= observed) / (1 + valid permutations).
  * Sequential stopping: after 99 permutations, only the genes with a permutation p-value
    (group effect or any contrast) below 0.1 continue to `n_permutations`. Stopped genes keep
    their 99-permutation p-values (resolution 0.01, all >= 0.1), so decisions at FDR 0.05 are
    not affected.
  * When there are no more than `n_permutations` distinct arrangements of the labels, all of them
    are used (exact test) and a warning is given.
  * `iDESC()` stops with an error if the permutation test is to be used and a covariate varies
    within subject or a subject has cells in several groups.
* **One per-gene status.** The `Status` column of `model_results` (added after 0.2.0) and the
  new `gene_status` element use the same values: `ok`, `no_standard_errors` (fitted, but a group
  coefficient has no Wald p-value), `fit_error`, `one_group` (expressed in fewer than two
  groups) and `filtered` (removed by the gene filter; only in `gene_status`). `gene_status` has
  one row per input gene, with the optimizer `convergence` code, `pdHess`, `convergence_null`,
  `n_perm_failed` and the error `message`. `model_results` has one row per gene that passed the
  gene filter (genes that were not fitted have `NA`), and
  `problematic_genes_cells_expressed$Reason` is the same status. The earlier values `OK`,
  `model_crash`, `convergence_error`, `invalid_likelihood` and `invalid_pval` are not used any
  more.
* `fit_stats = TRUE` adds `LogLik`, `Npar` and `AIC`.
* New element `permutation` in the result (whether the test was used, subjects per group,
  number of permutations, genes that were not stopped, ...).
* Unit tests (testthat).

## Changed

* Authors: Yunqing Liu (original author), Pablo Villoslada-Blanco (author and maintainer of this
  fork) and Adrián Santiso (contributor).
* The last column is `Minus2LogLik` (`Deviance` in 0.2.0, which was always `NA` with recent
  versions of glmmTMB); it is computed from the optimizer objective and is therefore available
  for every fitted gene.
* The coefficient columns keep the names of 0.2.0, taken from the model (`Beta_group<level>`,
  `Beta_<covariate>`); the unreleased commits after 0.2.0 had switched to
  `Beta_<group_var><level>`, assigned by position.
* FDR columns are BH adjustments across all genes with a p-value (the unreleased commits
  restricted them to genes with status `OK`).
* A full model that fails, or that ends with a lower likelihood than its null model, is fitted
  once more from the null estimates and the better fit is kept.
* `norm_opt` is matched with `match.arg()`, so the documented default (`"SeqDepth"`) works when
  the argument is not given, and a wrong value is an error.
* `iDESC()` checks its input: variables present in `meta`, no missing values in group, subject
  or covariates, positive normalization factors.
* The group is always used as a factor without unused levels; character covariates are
  converted to factors.

## Bug fixes

* `preprocessing()`: the proportion of expressing cells per subject was computed with the
  positions of the filtered metadata applied to the unfiltered barcodes, i.e. on the wrong
  cells whenever a cell had been removed. The gene filter (`gene_sub_filtering > 0`) now works
  as documented: a gene is kept if, in at least one group, it is expressed in at least
  `gene_cell_filtering` of the cells in at least `gene_sub_filtering` x (number of subjects of
  the group) subjects. The branch for groups with more than three levels no longer fails when
  a subject has been removed, and groups with a single subject are supported. `iDESC()` passes
  `cores` to `preprocessing()`.
* Genes whose model failed were collected inside a forked `mclapply()` worker and never reached
  the main process. Failures are now returned by the worker and reported (`fit_error`), the
  extraction of the coefficients is protected by `tryCatch()`, and `try-error` elements returned
  by `mclapply()` are handled, so one bad gene cannot stop or corrupt the run.
* Results of a gene were assigned to the output columns by position; they are now matched by
  coefficient name.
