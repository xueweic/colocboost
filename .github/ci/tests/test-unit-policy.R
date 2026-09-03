#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(testthat))

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
runner <- file.path(root, ".github", "ci", "run-unit-tests.R")
fixture_package <- file.path(
  root, ".github", "ci", "tests", "fixtures", "unit", "package"
)
rscript <- file.path(R.home("bin"), "Rscript")

make_waiver <- function(
  id = "fixture-waiver",
  context = "r-cmd-check-installed",
  file = "test-allowed-skip.R",
  test_title = "allowed skip case",
  reason = "allowed fixture reason",
  expected_count = 1L,
  expires = "2099-12-31"
) {
  list(
    id = id,
    context = context,
    file = file,
    test_title = test_title,
    reason = reason,
    expected_count = expected_count,
    rationale = "Synthetic policy fixture.",
    expires = expires
  )
}

write_policy <- function(waivers = list()) {
  path <- tempfile(fileext = ".yml")
  yaml::write_yaml(
    list(version = 1L, unit_tests = list(allowed_skips = waivers)),
    path
  )
  path
}

empty_policy <- write_policy()

install_library <- tempfile("unit-policy-library-")
dir.create(install_library)
install_output <- system2(
  file.path(R.home("bin"), "R"),
  c(
    "CMD", "INSTALL", "--install-tests", "--no-byte-compile", "--no-staged-install",
    paste0("--library=", shQuote(install_library)), shQuote(fixture_package)
  ),
  stdout = TRUE,
  stderr = TRUE
)
install_status <- attr(install_output, "status")
if (is.null(install_status)) {
  install_status <- 0L
}
if (install_status != 0L) {
  stop(paste(c("Unable to install synthetic fixture package:", install_output), collapse = "\n"))
}

run_runner <- function(
  filter = NULL,
  context = "local-r44",
  policy = empty_policy,
  load_package = "source",
  extra_args = character(),
  include_defaults = TRUE
) {
  output <- tempfile(fileext = ".json")
  default_args <- c(
    paste0("--package=", fixture_package),
    paste0("--load-package=", load_package),
    paste0("--context=", context),
    paste0("--policy=", policy),
    paste0("--output=", output)
  )
  if (!is.null(filter)) {
    default_args <- c(default_args, paste0("--filter=", filter))
  }
  args <- c(runner, if (include_defaults) default_args, extra_args)
  command_output <- suppressWarnings(system2(
    rscript,
    args,
    stdout = TRUE,
    stderr = TRUE,
    env = paste0("R_LIBS=", install_library)
  ))
  status <- attr(command_output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  list(
    status = unname(status),
    output = command_output,
    path = output,
    report = if (file.exists(output)) {
      jsonlite::fromJSON(output, simplifyVector = FALSE)
    } else {
      NULL
    }
  )
}

expect_summary <- function(report, ...) {
  expected <- list(...)
  for (name in names(expected)) {
    expect_identical(report$summary[[name]], expected[[name]], info = name)
  }
}

test_that("a passing fixture produces a structured passing report", {
  result <- run_runner("pass")

  expect_identical(result$status, 0L)
  expect_identical(result$report$schema_version, 1L)
  expect_identical(result$report$kind, "unit-tests")
  expect_identical(result$report$context, "local-r44")
  expect_identical(result$report$status, "pass")
  expect_summary(
    result$report,
    total = 1L, pass = 1L, failure = 0L, error = 0L, warning = 0L,
    empty = 0L, skip_like_success = 0L, skip = 0L
  )
  expect_identical(result$report$cases[[1]]$class, "pass")
  expect_identical(result$report$cases[[1]]$file, "test-pass.R")
  expect_identical(result$report$cases[[1]]$test_title, "pass case")
  expect_true("message" %in% names(result$report$cases[[1]]))
  expect_true("call" %in% names(result$report$cases[[1]]))
  expect_true("line" %in% names(result$report$cases[[1]]))
  expect_length(result$report$violations, 0L)
})

test_that("failure, error, and warning each fail with their strict class", {
  expectations <- list(
    failure = list(field = "failure", message = "Expected 1 to equal 2"),
    error = list(field = "error", message = "fixture error"),
    warning = list(field = "warning", message = "fixture warning")
  )

  for (filter in names(expectations)) {
    result <- run_runner(filter)
    expected <- expectations[[filter]]
    expect_identical(result$status, 1L, info = filter)
    expect_identical(result$report$status, "fail", info = filter)
    expect_identical(result$report$summary[[expected$field]], 1L, info = filter)
    expect_identical(result$report$cases[[1]]$class, filter, info = filter)
    expect_match(result$report$cases[[1]]$message, expected$message, info = filter)
  }
})

test_that("empty and executed skip-like success are never accepted", {
  empty <- run_runner("empty")
  skip_like <- run_runner("skip-like")

  expect_identical(empty$status, 1L)
  expect_summary(empty$report, empty = 1L, skip = 0L)
  expect_identical(empty$report$cases[[1]]$class, "empty")
  expect_identical(skip_like$status, 1L)
  expect_summary(skip_like$report, skip_like_success = 1L, pass = 0L)
  expect_identical(skip_like$report$cases[[1]]$class, "skip-like-success")
})

test_that("source mode rejects skips even if a source-context waiver is present", {
  policy <- write_policy(list(make_waiver(context = "local-r44")))
  result <- run_runner("allowed-skip", policy = policy)

  expect_identical(result$status, 1L)
  expect_summary(result$report, skip = 1L)
  expect_length(result$report$matched_allowances, 0L)
  expect_length(result$report$unused_allowances, 0L)
  expect_true(any(vapply(
    result$report$violations,
    function(x) identical(x$type, "unexpected-skip"),
    logical(1)
  )))
})

test_that("an exact installed-context waiver accepts exactly one skip", {
  policy <- write_policy(list(make_waiver()))
  result <- run_runner(
    "allowed-skip",
    context = "r-cmd-check-installed",
    policy = policy,
    load_package = "installed"
  )

  expect_identical(result$status, 0L)
  expect_identical(result$report$status, "pass")
  expect_summary(result$report, total = 1L, skip = 1L)
  expect_length(result$report$matched_allowances, 1L)
  expect_identical(result$report$matched_allowances[[1]]$id, "fixture-waiver")
  expect_identical(result$report$matched_allowances[[1]]$observed_count, 1L)
  expect_length(result$report$unused_allowances, 0L)
})

test_that("unexpected and changed skip identity are rejected exactly", {
  unexpected <- run_runner(
    "unexpected-skip",
    context = "r-cmd-check-installed",
    policy = write_policy(),
    load_package = "installed"
  )
  changed_reason <- run_runner(
    "allowed-skip",
    context = "r-cmd-check-installed",
    policy = write_policy(list(make_waiver(reason = "changed reason"))),
    load_package = "installed"
  )
  changed_context <- run_runner(
    "allowed-skip",
    context = "r-cmd-check-installed",
    policy = write_policy(list(make_waiver(context = "other-installed-context"))),
    load_package = "installed"
  )

  expect_identical(unexpected$status, 1L)
  expect_identical(changed_reason$status, 1L)
  expect_length(changed_reason$report$unused_allowances, 1L)
  expect_identical(changed_context$status, 1L)
  expect_length(changed_context$report$unused_allowances, 0L)
  expect_true(all(vapply(
    list(unexpected, changed_reason, changed_context),
    function(x) any(vapply(
      x$report$violations,
      function(v) identical(v$type, "unexpected-skip"),
      logical(1)
    )),
    logical(1)
  )))
})

test_that("skip counts must match their allowance exactly", {
  policy <- write_policy(list(make_waiver(
    expected_count = 2L
  )))
  result <- run_runner(
    "allowed-skip",
    context = "r-cmd-check-installed",
    policy = policy,
    load_package = "installed"
  )

  expect_identical(result$status, 1L)
  expect_summary(result$report, total = 1L, skip = 1L)
  expect_true(any(vapply(
    result$report$violations,
    function(x) identical(x$type, "skip-count-mismatch"),
    logical(1)
  )))
})

test_that("unused and expired same-context waivers fail", {
  unused <- run_runner(
    "pass",
    context = "r-cmd-check-installed",
    policy = write_policy(list(make_waiver())),
    load_package = "installed"
  )
  expired <- run_runner(
    "allowed-skip",
    context = "r-cmd-check-installed",
    policy = write_policy(list(make_waiver(expires = "2000-01-01"))),
    load_package = "installed"
  )

  expect_identical(unused$status, 1L)
  expect_length(unused$report$unused_allowances, 1L)
  expect_true(any(vapply(
    unused$report$violations,
    function(x) identical(x$type, "unused-allowance"),
    logical(1)
  )))
  expect_identical(expired$status, 1L)
  expect_true(any(vapply(
    expired$report$violations,
    function(x) identical(x$type, "expired-allowance"),
    logical(1)
  )))
})

test_that("duplicate waivers are a policy infrastructure error", {
  waiver <- make_waiver()
  result <- run_runner(
    "allowed-skip",
    context = "r-cmd-check-installed",
    policy = write_policy(list(waiver, waiver)),
    load_package = "installed"
  )

  expect_identical(result$status, 2L)
  expect_identical(result$report$status, "error")
  expect_true(any(vapply(
    result$report$violations,
    function(x) identical(x$type, "policy-error"),
    logical(1)
  )))
})

test_that("malformed CLI variants exit 2 and write diagnostics when possible", {
  output <- tempfile(fileext = ".json")
  base <- c(
    paste0("--package=", fixture_package),
    "--load-package=source",
    "--context=local-r44",
    paste0("--policy=", empty_policy),
    paste0("--output=", output)
  )
  variants <- list(
    unknown = c(base, "--unknown=value"),
    duplicate = c(base, "--context=duplicate"),
    malformed = c(base, "--filter"),
    missing = base[!grepl("^--policy=", base)]
  )

  for (name in names(variants)) {
    unlink(output)
    command_output <- suppressWarnings(system2(
      rscript,
      c(runner, variants[[name]]),
      stdout = TRUE,
      stderr = TRUE,
      env = paste0("R_LIBS=", install_library)
    ))
    status <- attr(command_output, "status")
    if (is.null(status)) status <- 0L
    expect_identical(unname(status), 2L, info = name)
    expect_true(file.exists(output), info = name)
    report <- jsonlite::fromJSON(output, simplifyVector = FALSE)
    expect_identical(report$status, "error", info = name)
  }
})

test_that("exit-1 observation paths preserve valid diagnostic JSON", {
  result <- run_runner("failure")

  expect_identical(result$status, 1L)
  expect_true(file.exists(result$path))
  expect_silent(jsonlite::validate(result$path))
  expect_identical(result$report$status, "fail")
  expect_true(length(result$report$violations) > 0L)
})
