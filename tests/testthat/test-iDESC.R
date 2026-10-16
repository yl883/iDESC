res_default <- run_idesc(sim3)

direct_fit <- function(dat, gene, zi, full = TRUE) {
  d <- data.frame(y = dat$mat[gene, ], sf = 1, group = dat$meta$group, sub = factor(dat$meta$subject))
  form <- if (full) y ~ group + offset(log(sf)) + (1 | sub) else y ~ 1 + offset(log(sf)) + (1 | sub)
  suppressWarnings(glmmTMB::glmmTMB(form, data = d, family = glmmTMB::nbinom2, ziformula = zi))
}

test_that("the result has the documented structure", {
  expect_identical(names(res_default), c("model_results", "problematic_genes_cells_expressed",
                                         "total_cells", "gene_status", "permutation"))
  mr <- res_default$model_results
  expect_identical(colnames(mr), cols_030(c("groupB", "groupC"), pairs = "Beta_C_vs_B"))
  expect_identical(rownames(mr), rownames(sim3$mat))
  expect_true(all(mr$Status == "ok"))
  expect_null(res_default$problematic_genes_cells_expressed)
  expect_equal(unlist(res_default$total_cells), c(A = 120, B = 120, C = 120))
  expect_false(res_default$permutation$used)
  expect_equal(mr$FDR_Beta_groupC, stats::p.adjust(mr$Pval_Beta_groupC, method = "BH"))
})

test_that("covariates add their coefficients after the group", {
  res <- run_idesc(sim3, covariates = "x")
  expect_identical(colnames(res$model_results),
                   cols_030(c("groupB", "groupC", "x"), pairs = "Beta_C_vs_B", group_coefs = c("groupB", "groupC")))
  dat <- sim3
  dat$meta$batch <- ifelse(as.integer(factor(dat$meta$subject)) %% 2 == 0, "b2", "b1")
  res2 <- run_idesc(dat, covariates = c("x", "batch"), pairwise = FALSE)
  expect_true(all(c("Beta_x", "Pval_Beta_x", "FDR_Beta_x", "Beta_batchb2") %in% colnames(res2$model_results)))
  expect_error(run_idesc(sim3, covariates = "nope"), "Not found in `meta`")
})

test_that("zi_model = 'none' (default) is the negative binomial mixed model", {
  mr <- res_default$model_results
  expect_true(all(is.na(mr$Theta)))
  expect_identical(run_idesc(sim3, zi_model = "none")$model_results, mr)
  f <- direct_fit(sim3, "g01", ~0)
  expect_equal(mr["g01", "Beta_groupC"], unname(glmmTMB::fixef(f)$cond["groupC"]), tolerance = 1e-6)
  expect_equal(mr["g01", "Alpha"], unname(glmmTMB::fixef(f)$cond["(Intercept)"]), tolerance = 1e-6)
  expect_equal(mr["g01", "Minus2LogLik"], 2 * f$fit$objective, tolerance = 1e-8)
  expect_equal(mr["g01", "Dispersion"], 1 / stats::sigma(f), tolerance = 1e-5)
  expect_equal(mr["g01", "Sigma2"], unname(unlist(glmmTMB::VarCorr(f)$cond$sub)[1]), tolerance = 1e-5)
})

test_that("zi_model = 'constant' fits one zero-inflation intercept", {
  res <- run_idesc(sim3, zi_model = "constant")
  mr <- res$model_results
  expect_identical(colnames(mr), colnames(res_default$model_results))
  expect_equal(nrow(mr), nrow(sim3$mat))
  expect_true(all(is.finite(mr$Theta)))
  f <- direct_fit(sim3, "g01", ~1)
  expect_equal(mr["g01", "Beta_groupC"], unname(glmmTMB::fixef(f)$cond["groupC"]), tolerance = 1e-5)
  expect_equal(mr["g01", "Theta"], unname(glmmTMB::fixef(f)$zi[1]), tolerance = 1e-4)
})

test_that("the original zero-inflation model and its arguments are gone", {
  expect_error(run_idesc(sim3, zi_model = "idesc"), "should be one of")
  expect_error(run_idesc(sim3, span = 0.5), "unused argument")
  expect_error(run_idesc(sim3, loess_control = "control"), "unused argument")
  expect_error(run_idesc(sim3, use_ind_zero = TRUE), "unused argument")
  expect_false("zp_prediction" %in% getNamespaceExports("iDESC"))
  expect_false(any(c("Sigma2_ZI", "Deviance") %in% colnames(res_default$model_results)))
})

test_that("norm_opt has a working default and is checked", {
  res <- suppressWarnings(suppressMessages(
    iDESC(sim3$mat[1:3, ], sim3$meta, subject_var = "subject", group_var = "group", permutation = "no")
  ))
  res2 <- suppressWarnings(suppressMessages(
    iDESC(sim3$mat[1:3, ], sim3$meta, subject_var = "subject", group_var = "group", norm_opt = "SeqDepth",
          permutation = "no")
  ))
  expect_identical(res, res2)
  expect_error(run_idesc(sim3, norm_opt = "nope"))
})

test_that("the likelihood-ratio test of the group effect is always reported", {
  for (zi in c("none", "constant")) {
    mr <- run_idesc(sim3, zi_model = zi)$model_results
    expect_true(all(is.finite(mr$LRT_stat)))
    expect_true(all(mr$LRT_stat >= 0))
    expect_equal(mr$Pval_LRT_group, stats::pchisq(mr$LRT_stat, df = 2, lower.tail = FALSE))
    expect_equal(mr$FDR_LRT_group, stats::p.adjust(mr$Pval_LRT_group, method = "BH"))
    expect_lt(mr["g01", "Pval_LRT_group"], 1e-4)
    expect_lt(mr["g01", "FDR_LRT_group"], 1e-3)
    zif <- if (zi == "none") ~0 else ~1
    f1 <- direct_fit(sim3, "g01", zif)
    f0 <- direct_fit(sim3, "g01", zif, full = FALSE)
    expect_equal(mr["g01", "LRT_stat"], 2 * (f0$fit$objective - f1$fit$objective), tolerance = 1e-5)
  }
  # two groups: one degree of freedom
  dat2 <- sim_counts(n_groups = 2, n_genes = 4, seed = 2)
  mr2 <- run_idesc(dat2)$model_results
  expect_equal(mr2$Pval_LRT_group, stats::pchisq(mr2$LRT_stat, df = 1, lower.tail = FALSE))
})

test_that("pairwise contrasts between non-reference levels", {
  mr <- res_default$model_results
  expect_equal(mr$Beta_C_vs_B, mr$Beta_groupC - mr$Beta_groupB)
  expect_equal(mr$FDR_Beta_C_vs_B, stats::p.adjust(mr$Pval_Beta_C_vs_B, method = "BH"))
  expect_lt(mr["g01", "Pval_Beta_C_vs_B"], 1e-4)
  # same test as refitting with B as the reference level
  dat2 <- sim3
  dat2$meta$group <- stats::relevel(dat2$meta$group, ref = "B")
  res2 <- run_idesc(dat2)
  expect_equal(mr$Beta_C_vs_B, res2$model_results$Beta_groupC, tolerance = 1e-3)
  expect_equal(log10(mr$Pval_Beta_C_vs_B), log10(res2$model_results$Pval_Beta_groupC), tolerance = 0.02)
  # switched off
  mr_off <- run_idesc(sim3, pairwise = FALSE)$model_results
  expect_identical(colnames(mr_off), cols_030(c("groupB", "groupC")))
  expect_equal(mr_off, mr[, colnames(mr_off)])
})

test_that("four groups give three pairwise contrasts; two groups none; names are syntactic", {
  dat4 <- sim_counts(n_groups = 4, subj_per_group = 2, n_genes = 3, seed = 4)
  mr4 <- run_idesc(dat4)$model_results
  expect_identical(colnames(mr4), cols_030(c("groupB", "groupC", "groupD"),
                                           pairs = c("Beta_C_vs_B", "Beta_D_vs_B", "Beta_D_vs_C")))
  expect_equal(mr4$Beta_D_vs_C, mr4$Beta_groupD - mr4$Beta_groupC)

  dat2 <- sim_counts(n_groups = 2, n_genes = 4, seed = 2)
  expect_identical(colnames(run_idesc(dat2, pairwise = TRUE)$model_results), cols_030("groupB"))

  dat3 <- sim3
  dat3$mat <- dat3$mat[c(1, 14), ]
  levels(dat3$meta$group) <- c("ctrl", "type 1", "type-2")
  mr3 <- run_idesc(dat3)$model_results
  expect_true(all(c("Beta_grouptype.1", "Beta_grouptype.2", "Beta_type.2_vs_type.1",
                    "Pval_Beta_type.2_vs_type.1", "FDR_Beta_type.2_vs_type.1") %in% colnames(mr3)))
  expect_equal(mr3$Beta_type.2_vs_type.1, mr3$Beta_grouptype.2 - mr3$Beta_grouptype.1)
})

test_that("a group that is not a factor is used as a factor", {
  dat <- sim3
  dat$meta$group <- as.character(dat$meta$group)
  expect_equal(run_idesc(dat), res_default)
})

test_that("fit_stats = TRUE adds LogLik, Npar and AIC", {
  mr <- run_idesc(sim3, zi_model = "constant", fit_stats = TRUE)$model_results
  expect_identical(tail(colnames(mr), 3), c("LogLik", "Npar", "AIC"))
  expect_true(all(mr$Npar == 6))
  expect_equal(mr$AIC, -2 * mr$LogLik + 2 * mr$Npar)
  expect_equal(mr$Minus2LogLik, -2 * mr$LogLik)
})

test_that("gene_status has one row per input gene and uses the statuses of model_results", {
  gs <- res_default$gene_status
  expect_identical(colnames(gs), c("gene", "status", "zi_status", "convergence", "pdHess", "convergence_null",
                                   "n_perm_failed", "message"))
  expect_identical(gs$gene, rownames(sim3$mat))
  expect_identical(gs$status, res_default$model_results$Status)
  expect_true(all(gs$convergence %in% 0:1))
  expect_true(all(gs$convergence_null %in% 0:1))
  expect_true(is.logical(gs$pdHess) && !any(is.na(gs$pdHess)))
  expect_true(all(is.na(gs$n_perm_failed)))

  dat <- sim3
  # a gene expressed in one group only
  dat$mat["g03", dat$meta$group != "A"] <- 0L
  dat$mat["g03", which(dat$meta$group == "A")[1:5]] <- 1L
  # a gene that the gene filter removes: expressed in a single cell
  dat$mat["g05", ] <- 0L
  dat$mat["g05", 1] <- 1L
  # a gene whose model cannot be fitted (negative values are rejected by the family)
  dat$mat["g07", which(dat$mat["g07", ] > 0)[1]] <- -1L
  res <- run_idesc(dat, gene_sub_filtering = 0.5, gene_cell_filtering = 0.05, cores = 2)
  gs <- res$gene_status
  mr <- res$model_results
  expect_identical(gs$gene, rownames(dat$mat))
  st <- stats::setNames(gs$status, gs$gene)
  expect_identical(unname(st[c("g03", "g05", "g07")]), c("one_group", "filtered", "fit_error"))
  expect_true(is.character(gs$message[gs$gene == "g07"]) && nzchar(gs$message[gs$gene == "g07"]))
  expect_true(all(is.na(gs$message[gs$gene != "g07"])))
  expect_true(all(st[setdiff(names(st), c("g03", "g05", "g07"))] %in% c("ok", "no_standard_errors")))
  expect_true(all(is.na(gs$convergence[!gs$status %in% c("ok", "no_standard_errors")])))
  expect_true(all(is.na(gs$pdHess[!gs$status %in% c("ok", "no_standard_errors")])))
  # model_results: every gene that passed the gene filter, with the same status
  expect_identical(rownames(mr), setdiff(rownames(dat$mat), "g05"))
  expect_identical(mr$Status, unname(st[rownames(mr)]))
  expect_true(all(is.na(mr[c("g03", "g07"), -1])))
  expect_true(all(is.finite(mr[setdiff(rownames(mr), c("g03", "g07")), "LRT_stat"])))
  # FDR over the genes with a p-value only
  ok <- !is.na(mr$Pval_LRT_group)
  expect_equal(sum(ok), nrow(mr) - 2)
  expect_equal(mr$FDR_LRT_group[ok], stats::p.adjust(mr$Pval_LRT_group[ok], method = "BH"))
  # the failed gene does not break the run, and is reported (bug fix)
  prob <- res$problematic_genes_cells_expressed
  expect_identical(colnames(prob), c("A", "B", "C", "Reason"))
  expect_identical(prob[c("g03", "g07"), "Reason"], c("one_group", "fit_error"))
  expect_identical(rownames(prob), rownames(mr)[mr$Status != "ok"])
  expect_equal(unlist(prob["g03", c("A", "B", "C")]),
               c(A = sum(dat$mat["g03", dat$meta$group == "A"] > 0), B = 0, C = 0))
})

test_that("a sparse count matrix gives the same result as a dense one", {
  sp <- list(mat = Matrix::Matrix(sim3$mat, sparse = TRUE), meta = sim3$meta)
  expect_s4_class(sp$mat, "dgCMatrix")
  expect_equal(run_idesc(sp), res_default)
})
