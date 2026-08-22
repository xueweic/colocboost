library(testthat)

test_that("colocalization candidate grouping supports more than 10000-byte signatures", {
  n_outcomes <- 4589L
  update <- matrix(-1L, nrow = n_outcomes, ncol = 3L)
  update[seq_len(100L), 1:2] <- 1L
  update[101:200, 3] <- 1L

  old_signature <- paste0(update[, 1], collapse = ",")
  expect_gt(nchar(old_signature, type = "bytes"), 10000L)

  result <- .group_coloc_candidates(update, seq_len(ncol(update)))

  expect_equal(unname(result), list(c(1L, 2L), 3L))
  expect_true(all(nchar(names(result), type = "bytes") < 10000L))
})
