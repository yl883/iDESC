# Small simulated sparse count data: mostly zeros, subjects nested in groups.
# Gene "g01" is strongly differentially expressed in the last group.
sim_counts <- function(n_groups = 3, subj_per_group = 3, cells_per_subj = 40, n_genes = 14, seed = 1) {
  set.seed(seed)
  grp_names <- LETTERS[seq_len(n_groups)]
  n_subj <- n_groups * subj_per_group
  subject <- rep(sprintf("s%02d", seq_len(n_subj)), each = cells_per_subj)
  group <- rep(rep(grp_names, each = subj_per_group), each = cells_per_subj)
  n_cells <- length(subject)
  base_mu <- seq(0.08, 0.6, length.out = n_genes)
  sub_eff <- matrix(stats::rnorm(n_genes * n_subj, sd = 0.5), n_genes, n_subj)
  sub_idx <- as.integer(factor(subject))
  mat <- t(sapply(seq_len(n_genes), function(g) {
    mu <- base_mu[g] * exp(sub_eff[g, sub_idx])
    if (g == 1) mu <- ifelse(group == grp_names[n_groups], mu * 25, mu)
    stats::rnbinom(n_cells, mu = mu, size = 2)
  }))
  # every cell expresses at least one gene, so that no cell is dropped
  mat[n_genes, colSums(mat) == 0] <- 1L
  rownames(mat) <- sprintf("g%02d", seq_len(n_genes))
  colnames(mat) <- sprintf("c%04d", seq_len(n_cells))
  meta <- data.frame(subject = subject, group = factor(group, levels = grp_names),
                     x = rep(stats::rnorm(n_subj), each = cells_per_subj),
                     row.names = colnames(mat), stringsAsFactors = FALSE)
  list(mat = mat, meta = meta)
}

# iDESC() without messages; the permutation test is off unless asked for
run_idesc <- function(dat, ..., permutation = "no") {
  suppressWarnings(suppressMessages(
    iDESC(dat$mat, dat$meta, subject_var = "subject", group_var = "group", norm_opt = "None",
          permutation = permutation, ...)
  ))
}

# Expected columns of model_results
cols_030 <- function(coefs, pairs = NULL, perm = FALSE, fit_stats = FALSE, group_coefs = coefs) {
  triple <- function(x) if (length(x) == 0) character(0) else as.vector(rbind(x, paste0("Pval_", x), paste0("FDR_", x)))
  contrasts <- c(paste0("Beta_", group_coefs), pairs)
  c("Status", "zi_status", "Alpha", triple(paste0("Beta_", coefs)), triple(pairs),
    "LRT_stat", "Pval_LRT_group", "FDR_LRT_group",
    if (perm) c("Pval_perm_group", "FDR_perm_group",
                as.vector(rbind(paste0("Pval_perm_", contrasts), paste0("FDR_perm_", contrasts))), "N_perm"),
    "Dispersion", "Sigma2", "Theta", "Minus2LogLik",
    if (fit_stats) c("LogLik", "Npar", "AIC"))
}

sim3 <- sim_counts()
# three genes only (g01 is the differentially expressed one), for the permutation tests
sim_small <- sim_counts(cells_per_subj = 30, n_genes = 3, seed = 5)
