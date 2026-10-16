# Two groups of four subjects (10 cells each) plus three small subjects that are removed.
# Cells are in random order, so that positions in the filtered meta differ from the original ones.
make_filter_data <- function() {
  set.seed(42)
  subject <- c(rep(paste0("small", 1:3), each = 4), rep(paste0("s", 1:8), each = 10))
  group <- ifelse(subject %in% c("small1", "small2", paste0("s", 1:4)), "A", "B")
  n <- length(subject)
  cell_in_sub <- stats::ave(seq_len(n), subject, FUN = seq_along)
  genes <- c("base", "keepA", "keepB", "drop_one_subject", "drop_few_cells", "drop_small_only")
  mat <- matrix(0L, length(genes), n, dimnames = list(genes, paste0("cell", seq_len(n))))
  mat["base", ] <- 1L                                                   # all cells
  mat["keepA", subject %in% c("s1", "s2")] <- 3L                        # 2 of 4 subjects of A
  mat["keepB", subject %in% c("s5", "s6", "s7") & cell_in_sub <= 6] <- 2L  # 60% of cells, 3 of 4 subjects of B
  mat["drop_one_subject", subject %in% c("s1", "s5")] <- 1L             # 1 of 4 subjects in each group
  mat["drop_few_cells", cell_in_sub <= 4 & !grepl("small", subject)] <- 1L  # 40% of cells in every subject
  mat["drop_small_only", grepl("small", subject)] <- 5L                 # only in removed subjects
  meta <- data.frame(subject = subject, group = group, row.names = colnames(mat), stringsAsFactors = FALSE)
  ord <- sample(n)
  list(mat = mat[, ord], meta = meta[ord, ])
}

test_that("gene filter keeps/drops the right genes when cells are removed", {
  d <- make_filter_data()
  out <- preprocessing(d$mat, d$meta, "subject", "group", sub_cell_filtering = 5,
                       gene_sub_filtering = 0.5, gene_cell_filtering = 0.5, ncell_filtering = 1)
  expect_identical(rownames(out$mat), c("base", "keepA", "keepB"))
  expect_equal(ncol(out$mat), 80)
  expect_identical(colnames(out$mat), rownames(out$meta))
  expect_false(any(grepl("small", out$meta$subject)))
  expect_identical(out$mat, d$mat[rownames(out$mat), colnames(out$mat)])
})

test_that("defaults of the gene filter keep every gene", {
  d <- make_filter_data()
  out <- preprocessing(d$mat, d$meta, "subject", "group", sub_cell_filtering = 5,
                       gene_sub_filtering = 0, gene_cell_filtering = 0.05, ncell_filtering = 1)
  expect_identical(rownames(out$mat), rownames(d$mat))
  expect_equal(ncol(out$mat), 80)
})

test_that("a group with more than 3 levels uses all subjects together", {
  d <- make_filter_data()
  d$meta$group <- paste0("G", (as.integer(factor(d$meta$subject)) %% 4) + 1)
  expect_output(
    out <- preprocessing(d$mat, d$meta, "subject", "group", sub_cell_filtering = 5,
                         gene_sub_filtering = 0.25, gene_cell_filtering = 0.5, ncell_filtering = 1),
    "more than 3 levels"
  )
  # 8 subjects remain: a gene needs 2 subjects with >= 50% of cells expressing it
  expect_identical(rownames(out$mat), c("base", "keepA", "keepB", "drop_one_subject"))
})

test_that("ncell_filtering removes cells with too few expressed genes", {
  d <- make_filter_data()
  out <- preprocessing(d$mat, d$meta, "subject", "group", sub_cell_filtering = 5,
                       gene_sub_filtering = 0.5, gene_cell_filtering = 0.5, ncell_filtering = 2)
  expect_true(all(colSums(out$mat != 0) >= 2))
  expect_lt(ncol(out$mat), 80)
})
