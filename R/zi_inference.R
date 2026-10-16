# Choose a fit for observed-data coefficient inference only. The caller retains the original
# ZINB fit for likelihood summaries, the group LRT and all permutation statistics.
zi_inference_fit <- function(original, fit, group_coefs, boundary_tol) {
  usable <- function(f) {
    tryCatch({
      ct <- summary(f)$coefficients$cond
      isTRUE(f$sdr$pdHess) && all(group_coefs %in% rownames(ct)) &&
        all(is.finite(ct[group_coefs, "Std. Error"])) &&
        all(is.finite(ct[group_coefs, "Pr(>|z|)"]))
    }, error = function(e) FALSE)
  }
  ll <- function(f) -as.numeric(f$fit$objective)[1]
  out <- list(fit = original, status = "unidentified")
  if (usable(original)) return(list(fit = original, status = "estimated"))
  nb <- tryCatch(fit(zi = ~ 0, se = TRUE), error = function(e) NULL)
  if (is.null(nb) || !is.finite(ll(original)) || !is.finite(ll(nb))) return(out)
  # Absolute difference avoids declaring a substantially better NB fit a boundary solution.
  # An NB fit without usable inference would discard any SEs of the original fit: go on to retries.
  if (2 * abs(ll(original) - ll(nb)) < boundary_tol && usable(nb)) {
    return(list(fit = nb, status = "boundary"))
  }
  start <- tryCatch({
    pl <- nb$obj$env$parList(nb$fit$par, nb$fit$parfull)
    c(list(beta = unname(glmmTMB::fixef(nb)$cond)),
      pl[intersect(c("betadisp", "betad"), names(pl))],
      list(theta = pl$theta, betazi = stats::qlogis(0.1)))
  }, error = function(e) NULL)
  if (is.null(start)) return(out)
  for (bfgs in c(FALSE, TRUE)) {
    control <- if (bfgs) glmmTMB::glmmTMBControl(optimizer = stats::optim,
                                               optArgs = list(method = "BFGS")) else NULL
    retry <- tryCatch(fit(se = TRUE, start = start, control = control), error = function(e) NULL)
    if (!is.null(retry) && usable(retry) && is.finite(ll(retry)) &&
        abs(ll(retry) - ll(original)) <= 1e-4 * max(1, abs(ll(original)))) {
      return(list(fit = retry, status = "refit"))
    }
  }
  out
}
