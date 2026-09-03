sentinel <- paste0(
  "One run returned no colocalization; ",
  "entropy-vs-uniform comparison skipped"
)

test_that("error outranks all lower classes", {
  succeed(sentinel)
  warning("mixed warning before error")
  expect_equal(1, 2)
  stop("mixed fixture error")
})

test_that("failure outranks warning and lower classes", {
  succeed(sentinel)
  warning("mixed warning before failure")
  expect_equal(1, 2)
  skip("empty test")
})

test_that("warning outranks empty and lower classes", {
  succeed(sentinel)
  warning("mixed warning before empty")
  skip("empty test")
})

test_that("empty outranks skip-like success", {
  succeed(sentinel)
  skip("empty test")
})

test_that("skip-like success outranks skip", {
  succeed(sentinel)
  skip("mixed ordinary skip")
})

test_that("skip outranks ordinary pass", {
  succeed("mixed ordinary pass before skip")
  skip("mixed ordinary skip")
})

test_that("pass remains pass", {
  succeed("mixed ordinary pass")
})
