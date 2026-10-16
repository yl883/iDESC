test_that("NB inference never invokes ZINB recovery", {
  testthat::local_mocked_bindings(zi_inference_fit = function(...) stop("unexpected recovery"))
  a <- run_idesc(sim_small, zi_model = "none", zi_boundary_tol = 0)
  b <- run_idesc(sim_small, zi_model = "none", zi_boundary_tol = 1)
  expect_identical(a, b)
  expect_true(all(a$model_results$zi_status == "none"))
  expect_true(all(a$gene_status$zi_status == "none"))
})

test_that("boundary inference uses NB but likelihood summaries retain ZINB", {
  dat <- sim3
  dat$mat <- rbind(dat$mat[c("g02", "g08"), , drop = FALSE], background = 1)
  # The background row keeps zero-count cells through preprocessing.
  nb <- run_idesc(dat, ncell_filtering = 0, zi_model = "none", fit_stats = TRUE)
  zi <- run_idesc(dat, ncell_filtering = 0, zi_model = "constant", fit_stats = TRUE)
  rows <- zi$model_results$zi_status == "boundary"
  expect_true(any(rows))
  expect_true(all(zi$model_results$Status[rows] == "ok"))
  cols <- c("Alpha", grep("^(Beta_|Pval_Beta_)", names(nb$model_results), value = TRUE), "Dispersion", "Sigma2")
  expect_equal(zi$model_results[rows, cols], nb$model_results[rows, cols], tolerance = 1e-4)
  testthat::local_mocked_bindings(zi_inference_fit = function(original, ...) list(fit = original, status = "unidentified"))
  original <- run_idesc(dat, ncell_filtering = 0, zi_model = "constant", fit_stats = TRUE)
  cols <- c("LRT_stat", "Pval_LRT_group", "LogLik", "Minus2LogLik", "Npar", "AIC", "Theta")
  expect_identical(zi$model_results[, cols], original$model_results[, cols])
  expect_identical(zi$gene_status$pdHess, original$gene_status$pdHess)
})

test_that("real zero inflation is estimated without fallback", {
  dat <- sim_counts(n_groups = 2, subj_per_group = 4, cells_per_subj = 80, n_genes = 2, seed = 41)
  set.seed(42)
  sub <- as.integer(factor(dat$meta$subject))
  mu <- 12 * exp(rnorm(8, sd = 0.5)[sub])
  dat$mat[1, ] <- rnbinom(length(mu), mu = mu, size = 4) * rbinom(length(mu), 1, 0.6)
  dat$mat[2, ] <- 1L
  res <- run_idesc(dat, zi_model = "constant", ncell_filtering = 0)$model_results[1, , drop = FALSE]
  expect_identical(res$zi_status, "estimated")
  expect_identical(res$Status, "ok")
  expect_equal(res$Theta, qlogis(0.4), tolerance = 0.25)
})

test_that("ZINB permutation statistics are reproducible and unaffected by recovery", {
  dat <- sim3
  dat$mat <- rbind(dat$mat[c("g02", "g08"), , drop = FALSE], background = 1)
  run <- function() run_idesc(dat, ncell_filtering = 0, zi_model = "constant", permutation = "yes",
                              n_permutations = 3, sequential = FALSE, seed = 71)
  a <- run()
  expect_true(any(a$model_results$zi_status == "boundary"))
  expect_identical(a, run())
  testthat::local_mocked_bindings(zi_inference_fit = function(original, ...) list(fit = original, status = "unidentified"))
  b <- run()
  cols <- grep("perm|LRT|Minus2LogLik", names(a$model_results), value = TRUE)
  expect_identical(a$model_results[, cols], b$model_results[, cols])
  expect_identical(a$gene_status$n_perm_failed, b$gene_status$n_perm_failed)
})

test_that("retries require usable inference and a matching original likelihood", {
  # Real fitted objects keep the extraction and starting-value paths under test; the fit
  # callback controls optimizer outcomes without relying on an unstable convergence failure.
  dat <- sim_small
  d <- data.frame(y = dat$mat[1, ], group = dat$meta$group, sub = factor(dat$meta$subject))
  good <- suppressWarnings(glmmTMB::glmmTMB(y ~ group + (1 | sub), data = d,
                                           family = glmmTMB::nbinom2, ziformula = ~0))
  expect_true(isTRUE(good$sdr$pdHess))
  original <- good
  original$sdr$pdHess <- FALSE
  nb <- good
  nb$fit$objective <- original$fit$objective + 1
  attempts <- list()
  callback <- function(...) {
    args <- list(...)
    attempts[[length(attempts) + 1L]] <<- args
    if (length(attempts) == 1L) return(nb)
    if (length(attempts) == 2L) return(original)
    good
  }
  r <- zi_inference_fit(original, callback, c("groupB", "groupC"), 1e-3)
  expect_identical(r$status, "refit")
  expect_identical(r$fit, good)
  expect_length(attempts, 3)
  expect_equal(attempts[[2]]$start$betazi, qlogis(0.1))
  expect_equal(attempts[[2]]$start$beta, unname(glmmTMB::fixef(nb)$cond))
  expect_true(any(c("betadisp", "betad") %in% names(attempts[[2]]$start)))
  expect_identical(attempts[[3]]$control$optimizer, stats::optim)
  expect_identical(attempts[[3]]$control$optArgs$method, "BFGS")
  for (delta in c(-1, 1)) {
    other <- good
    other$fit$objective <- original$fit$objective + delta
    r <- zi_inference_fit(original, function(...) other, c("groupB", "groupC"), 1e-3)
    expect_identical(r$status, "unidentified")
    expect_identical(r$fit, original)
  }
  # A boundary NB fit without a positive definite Hessian is not accepted; recovery goes on to retries
  same <- original
  calls <- 0L
  r <- zi_inference_fit(original, function(...) { calls <<- calls + 1L; same }, c("groupB", "groupC"), 1e-3)
  expect_identical(r$status, "unidentified")
  expect_identical(r$fit, original)
  expect_identical(calls, 3L)
  r <- zi_inference_fit(good, function(...) stop("should not refit"), c("groupB", "groupC"), 1e-3)
  expect_identical(r$status, "estimated")
  r <- zi_inference_fit(original, function(...) stop("fit failed"), c("groupB", "groupC"), 1e-3)
  expect_identical(r$status, "unidentified")
})
