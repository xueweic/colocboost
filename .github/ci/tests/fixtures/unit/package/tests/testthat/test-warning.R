test_that("warning case", {
  warning("fixture warning")
  succeed("warning fixture reached its final expectation")
})
