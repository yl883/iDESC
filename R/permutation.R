# Internal helpers of the permutation test of iDESC() (group labels permuted between subjects)

# Number of distinct ways of assigning the group labels to the subjects (multinomial coefficient)
n_arrangements <- function(counts) {
  round(exp(lgamma(sum(counts) + 1) - sum(lgamma(counts + 1))))
}

# All distinct arrangements of the labels 1..K with the given counts.
# Returns an integer matrix with one arrangement per column.
all_arrangements <- function(counts) {
  n <- sum(counts)
  total <- n_arrangements(counts)
  out <- matrix(0L, n, total)
  cur <- integer(n)
  k <- 0L
  rec <- function(pos, left) {
    if (pos > n) {
      k <<- k + 1L
      out[, k] <<- cur
      return(invisible(NULL))
    }
    for (lev in which(left > 0)) {
      cur[pos] <<- lev
      left2 <- left
      left2[lev] <- left2[lev] - 1
      rec(pos + 1L, left2)
    }
  }
  rec(1L, counts)
  out
}

# Run `expr` after set.seed(seed) and restore the state of the random number generator
with_seed <- function(seed, expr) {
  genv <- globalenv()
  if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
    old <- get(".Random.seed", envir = genv, inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = genv))
  } else {
    on.exit(if (exists(".Random.seed", envir = genv, inherits = FALSE)) rm(".Random.seed", envir = genv))
  }
  set.seed(seed)
  expr
}

# Permuted group labels of the subjects.
# `codes`: observed group (integer code) of every subject. Returns a list with `labels` (integer
# matrix, subjects x permutations), `exhaustive` and `n_arrangements`.
# If the number of distinct arrangements other than the observed one is not larger than
# `n_permutations`, all of them are used (in random order): the test is then exact.
permutation_labels <- function(codes, n_permutations, seed) {
  codes <- as.integer(codes)
  counts <- tabulate(codes)
  total <- n_arrangements(counts)
  if (total - 1 <= n_permutations) {
    arr <- all_arrangements(counts)
    arr <- arr[, colSums(arr != codes) > 0, drop = FALSE]
    ord <- with_seed(seed, sample.int(ncol(arr)))
    return(list(labels = arr[, ord, drop = FALSE], exhaustive = TRUE, n_arrangements = total))
  }
  n <- length(codes)
  labels <- with_seed(seed, vapply(seq_len(n_permutations), function(i) codes[sample.int(n)], integer(n)))
  list(labels = labels, exhaustive = FALSE, n_arrangements = total)
}

# Permutation p-values of one gene.
# `obs`: observed statistics; `stats`: matrix permutations x statistics (rows of permutations that
# were run). A permutation is valid when all statistics that are available for the observed data
# are finite. p = (1 + number of valid permutations with statistic >= observed) / (1 + valid).
permutation_pvalues <- function(obs, stats) {
  need <- is.finite(obs)
  valid <- if (nrow(stats) > 0) rowSums(!is.finite(stats[, need, drop = FALSE])) == 0 else logical(0)
  n_valid <- sum(valid)
  p <- rep(NA_real_, length(obs))
  if (n_valid > 0) {
    for (k in which(need)) {
      # small tolerance, so that an arrangement equivalent to the observed one counts as a tie
      tol <- 1e-8 + 1e-6 * abs(obs[k])
      p[k] <- (1 + sum(stats[valid, k] >= obs[k] - tol)) / (1 + n_valid)
    }
  }
  list(p = p, n_valid = n_valid, n_failed = sum(!valid))
}
