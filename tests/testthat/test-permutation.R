perm_cols <- function(mr) grep("perm", colnames(mr), value = TRUE)

res_perm <- run_idesc(sim_small, permutation = "yes", n_permutations = 49, cores = 1)

test_that("permutation columns are added and p-values are in (0, 1]", {
  mr <- res_perm$model_results
  expect_identical(colnames(mr), cols_030(c("groupB", "groupC"), pairs = "Beta_C_vs_B", perm = TRUE))
  pcols <- grep("^Pval_perm_", colnames(mr), value = TRUE)
  expect_identical(pcols, c("Pval_perm_group", "Pval_perm_Beta_groupB", "Pval_perm_Beta_groupC",
                            "Pval_perm_Beta_C_vs_B"))
  for (pc in pcols) {
    p <- mr[[pc]]
    expect_true(all(p > 0 & p <= 1))
    # (1 + b) / (1 + 49): multiples of 1/50
    expect_equal(p * 50, round(p * 50))
    expect_equal(mr[[sub("^Pval_", "FDR_", pc)]], stats::p.adjust(p, method = "BH"))
  }
  expect_identical(mr$N_perm, rep(49L, 3))
  expect_identical(res_perm$gene_status$n_perm_failed, rep(0L, 3))
  pm <- res_perm$permutation
  expect_true(pm$used)
  expect_false(pm$exhaustive)
  expect_equal(pm$n_permutations, 49)
  expect_equal(pm$n_arrangements, 1680)
  expect_equal(unname(pm$subjects_per_group), c(3, 3, 3))
})

test_that("the asymptotic part does not depend on the permutation test", {
  res_no <- run_idesc(sim_small)
  mr <- res_perm$model_results
  expect_equal(mr[, setdiff(colnames(mr), c(perm_cols(mr), "N_perm"))], res_no$model_results)
})

test_that("a strongly differentially expressed gene gets a small permutation p-value", {
  mr <- res_perm$model_results
  expect_lte(mr["g01", "Pval_perm_group"], 0.05)
  expect_lte(mr["g01", "Pval_perm_Beta_C_vs_B"], 0.1)
  expect_gt(min(mr[c("g02", "g03"), "Pval_perm_group"]), 0.05)
})

test_that("permutation results are reproducible and do not depend on cores", {
  expect_identical(run_idesc(sim_small, permutation = "yes", n_permutations = 49, cores = 1), res_perm)
  res2 <- run_idesc(sim_small, permutation = "yes", n_permutations = 49, cores = 2)
  expect_identical(res2$model_results[, perm_cols(res2$model_results)],
                   res_perm$model_results[, perm_cols(res_perm$model_results)])
  expect_equal(res2, res_perm)
  # another seed gives other permutations
  res3 <- run_idesc(sim_small, permutation = "yes", n_permutations = 49, cores = 2, seed = 2)
  mr3 <- res3$model_results
  expect_false(identical(mr3[, perm_cols(mr3)], res_perm$model_results[, perm_cols(mr3)]))
  expect_equal(mr3$LRT_stat, res_perm$model_results$LRT_stat)
})

test_that("the random number stream of the session is not touched", {
  set.seed(123)
  a <- stats::runif(1)
  set.seed(123)
  invisible(run_idesc(sim_small, permutation = "yes", n_permutations = 5, cores = 1))
  expect_identical(stats::runif(1), a)
})

test_that("sequential stopping: N_perm is 99 for stopped genes and n_permutations for the others", {
  res_seq <- run_idesc(sim_small, permutation = "yes", n_permutations = 119, cores = 2)
  res_all <- run_idesc(sim_small, permutation = "yes", n_permutations = 119, cores = 2, sequential = FALSE)
  mr <- res_seq$model_results
  pcols <- grep("^Pval_perm_", colnames(mr), value = TRUE)
  expect_true(all(mr$N_perm %in% c(99L, 119L)))
  expect_identical(unname(mr["g01", "N_perm"]), 119L)
  stopped <- mr$N_perm == 99L
  expect_true(any(stopped))
  expect_true(all(as.matrix(mr[stopped, pcols]) >= 0.1))
  expect_true(all(apply(as.matrix(mr[!stopped, pcols]), 1, min) < 0.1 | !stopped[!stopped]))
  expect_identical(res_seq$permutation$genes_all_permutations, rownames(mr)[!stopped])
  expect_equal(res_seq$permutation$first_block, 99)
  # without stopping every gene uses all permutations; genes that were not stopped are unchanged
  mr_all <- res_all$model_results
  expect_identical(mr_all$N_perm, rep(119L, 3))
  expect_identical(mr_all[!stopped, pcols], mr[!stopped, pcols])
  # a stopped gene would not have been significant anyway
  expect_true(all(as.matrix(mr_all[stopped, pcols]) > 0.05))
  # no stopping when n_permutations <= 99
  expect_identical(res_perm$permutation$genes_all_permutations, rownames(sim_small$mat))
})

test_that("permutation = 'auto' follows the number of subjects per group", {
  auto <- function(dat, ...) {
    suppressWarnings(suppressMessages(
      iDESC(dat$mat, dat$meta, subject_var = "subject", group_var = "group", norm_opt = "None",
            n_permutations = 9, cores = 2, ...)
    ))
  }
  dat <- sim_small
  dat$mat <- dat$mat[c(1, 3), ]
  # 3 subjects per group < 10: on
  r_on <- auto(dat)
  expect_true(r_on$permutation$used)
  expect_true("Pval_perm_group" %in% colnames(r_on$model_results))
  expect_identical(r_on$model_results$N_perm, rep(9L, 2))
  # threshold lowered to 3: off
  r_off <- auto(dat, min_subjects_asymptotic = 3)
  expect_false(r_off$permutation$used)
  expect_false(any(grepl("perm", colnames(r_off$model_results))))
  expect_false("N_perm" %in% colnames(r_off$model_results))
  # explicit choices override the rule
  expect_false(auto(dat, permutation = "no")$permutation$used)
  expect_true(auto(dat, permutation = "yes", min_subjects_asymptotic = 1)$permutation$used)
  # 10 subjects per group: off by default, on when the threshold is 11; unbalanced groups use the smallest
  big <- sim_counts(n_groups = 2, subj_per_group = 10, cells_per_subj = 12, n_genes = 2, seed = 7)
  expect_false(auto(big)$permutation$used)
  expect_true(auto(big, min_subjects_asymptotic = 11)$permutation$used)
  keep <- big$meta$subject != "s01"
  unbalanced <- list(mat = big$mat[, keep], meta = big$meta[keep, ])
  r_unb <- auto(unbalanced)
  expect_true(r_unb$permutation$used)
  expect_equal(unname(r_unb$permutation$subjects_per_group), c(9, 10))
  expect_message(
    suppressWarnings(iDESC(big$mat, big$meta, subject_var = "subject", group_var = "group",
                           norm_opt = "None", cores = 2)),
    "permutation p-values are not computed"
  )
})

test_that("permutation needs subject-level covariates and subjects nested in groups", {
  dat <- sim_small
  dat$mat <- dat$mat[c(1, 3), ]
  set.seed(11)
  dat$meta$cell_cov <- stats::rnorm(nrow(dat$meta))
  expect_error(run_idesc(dat, permutation = "yes", covariates = "cell_cov", n_permutations = 5),
               "`cell_cov` varies within subject")
  expect_error(run_idesc(dat, permutation = "auto", covariates = c("x", "cell_cov"), n_permutations = 5),
               "`cell_cov` varies within subject")
  # a subject-level covariate is fine, and the cell-level one is fine without permutation
  res <- run_idesc(dat, permutation = "yes", covariates = "x", n_permutations = 5, cores = 2)
  expect_identical(res$model_results$N_perm, rep(5L, 2))
  expect_true("Beta_cell_cov" %in% colnames(run_idesc(dat, permutation = "no", covariates = "cell_cov")$model_results))

  crossed <- dat
  first_cells <- !duplicated(crossed$meta$subject)
  crossed$meta$group[first_cells] <- ifelse(crossed$meta$group[first_cells] == "A", "B", "A")
  expect_error(run_idesc(crossed, permutation = "yes", n_permutations = 5), "several groups")
  expect_error(run_idesc(dat, permutation = "yes", n_permutations = 0), "n_permutations")
})

test_that("all arrangements are used when there are fewer than n_permutations", {
  dat <- sim_counts(n_groups = 2, subj_per_group = 2, cells_per_subj = 40, n_genes = 2, seed = 3)
  expect_warning(
    res <- suppressMessages(iDESC(dat$mat, dat$meta, subject_var = "subject", group_var = "group",
                                  norm_opt = "None", permutation = "yes", n_permutations = 49, cores = 2)),
    "Only 6 distinct arrangements"
  )
  expect_true(res$permutation$exhaustive)
  expect_equal(res$permutation$n_permutations, 5)
  expect_equal(res$permutation$n_permutations_requested, 49)
  expect_identical(res$model_results$N_perm, rep(5L, 2))
  # the arrangement with the two labels swapped gives the same statistics: p >= 2/6
  expect_true(all(res$model_results$Pval_perm_group >= 2 / 6 - 1e-12))
  expect_true(all(res$model_results$Pval_perm_Beta_groupB >= 2 / 6 - 1e-12))
})

test_that("permutation helpers", {
  arr <- iDESC:::all_arrangements(c(2, 1, 1))
  expect_equal(dim(arr), c(4, 12))
  expect_equal(ncol(unique(arr, MARGIN = 2)), 12)
  expect_true(all(apply(arr, 2, function(x) identical(tabulate(x, 3), c(2L, 1L, 1L)))))
  expect_equal(iDESC:::n_arrangements(c(5, 6, 6)), 5717712)

  codes <- rep(1:3, c(5, 6, 6))
  p1 <- iDESC:::permutation_labels(codes, 20, seed = 1)
  expect_identical(p1, iDESC:::permutation_labels(codes, 20, seed = 1))
  expect_false(identical(p1$labels, iDESC:::permutation_labels(codes, 20, seed = 2)$labels))
  expect_equal(dim(p1$labels), c(17, 20))
  expect_true(all(apply(p1$labels, 2, function(x) identical(tabulate(x, 3), c(5L, 6L, 6L)))))
  expect_false(p1$exhaustive)
  # the first permutations do not depend on how many are requested
  expect_identical(p1$labels[, 1:10], iDESC:::permutation_labels(codes, 10, seed = 1)$labels)

  ex <- iDESC:::permutation_labels(c(1, 1, 2, 2), 10, seed = 1)
  expect_true(ex$exhaustive)
  expect_equal(dim(ex$labels), c(4, 5))
  expect_false(any(apply(ex$labels, 2, function(x) identical(x, c(1L, 1L, 2L, 2L)))))

  # p = (1 + b) / (1 + valid); rows with a non-finite needed statistic are not valid
  stats <- cbind(c(1, 2, 3, NA, 5), c(0.1, 0.2, 0.3, 0.4, 0.5), NA)
  pv <- iDESC:::permutation_pvalues(c(3, 0.45, NA), stats)
  expect_equal(pv$n_valid, 4)
  expect_equal(pv$n_failed, 1)
  expect_equal(pv$p, c((1 + 2) / 5, (1 + 1) / 5, NA))
})
