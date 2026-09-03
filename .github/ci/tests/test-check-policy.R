#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(testthat))

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
runner <- file.path(root, ".github", "ci", "run-r-cmd-check.R")
rscript <- file.path(R.home("bin"), "Rscript")
python <- Sys.which("python")

source_sha <- paste(rep("a", 40L), collapse = "")
event_sha <- paste(rep("b", 40L), collapse = "")
tarball_sha256 <- paste(rep("c", 64L), collapse = "")

make_note_rule <- function(
  id = "installed-size-5mb",
  check = "checking installed package size",
  output = paste(
    "  installed size is  5.0Mb",
    "  sub-directories of 1Mb or more:",
    "    data   2.0Mb",
    "    doc    1.9Mb",
    sep = "\n"
  ),
  expires = "2099-12-31"
) {
  list(
    id = id,
    check = check,
    output = output,
    rationale = "Synthetic exact NOTE fixture.",
    expires = expires
  )
}

write_policy <- function(notes = list()) {
  path <- tempfile(fileext = ".yml")
  yaml::write_yaml(
    list(version = 1L, r_cmd_check = list(allowed_notes = notes)),
    path
  )
  path
}

empty_policy <- write_policy()

make_log <- function(chunks, status = "OK", done = TRUE, elapsed = FALSE) {
  lines <- c(
    "* using session charset: UTF-8",
    "* this is package 'colocboost' version '1.0.9'",
    chunks
  )
  if (done) {
    lines <- c(lines, "* DONE")
  }
  if (!is.null(status)) {
    lines <- c(lines, paste0("Status: ", status))
  }
  if (elapsed) {
    lines <- c(lines, "* elapsed time (check, wall clock): 0h 0m 1s")
  }
  paste(lines, collapse = "\n")
}

clean_chunks <- c(
  "* checking package dependencies ... OK",
  "This prose contains ERROR, WARNING, and NOTE but is not a heading status.",
  "* checking tests ... OK"
)

run_parser <- function(
  log_text,
  policy = empty_policy,
  check_exit_code = 0L,
  extra_args = character(),
  line_ending = "\n",
  trailing_newline = TRUE
) {
  log <- tempfile(fileext = ".log")
  output <- tempfile(fileext = ".json")
  serialized <- gsub("\n", line_ending, log_text, fixed = TRUE)
  if (trailing_newline) serialized <- paste0(serialized, line_ending)
  connection <- file(log, open = "wb")
  writeBin(charToRaw(serialized), connection)
  close(connection)
  args <- c(
    runner,
    paste0("--log=", log),
    "--reported-log-path=colocboost.Rcheck/00check.log",
    paste0("--policy=", policy),
    paste0("--output=", output),
    "--environment-id=fixture-env",
    paste0("--source-sha=", source_sha),
    paste0("--event-sha=", event_sha),
    paste0("--tarball-sha256=", tarball_sha256),
    paste0("--check-exit-code=", check_exit_code),
    extra_args
  )
  command_output <- suppressWarnings(system2(
    rscript,
    args,
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(command_output, "status")
  if (is.null(status)) status <- 0L
  list(
    status = unname(status),
    command_output = command_output,
    path = output,
    document = if (file.exists(output)) {
      jsonlite::fromJSON(output, simplifyVector = FALSE)
    } else {
      NULL
    }
  )
}

expect_contract_valid <- function(path) {
  code <- paste(
    "import json,sys",
    "sys.path.insert(0,'.github/ci')",
    "from result_contract import validate_result",
    "validate_result(json.load(open(sys.argv[1], encoding='utf-8')))",
    sep = ";"
  )
  validation <- system2(
    python,
    c("-B", "-c", shQuote(code), shQuote(path)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(validation, "status")
  if (is.null(status)) status <- 0L
  expect_identical(unname(status), 0L, info = paste(validation, collapse = "\n"))
}

expect_infrastructure_failure <- function(result) {
  expect_identical(result$status, 2L)
  expect_true(file.exists(result$path))
  expect_contract_valid(result$path)
  expect_identical(result$document$result_kind, "infrastructure")
  expect_identical(result$document$status, "fail")
}

test_that("clean headings and prose keywords produce a package pass", {
  result <- run_parser(make_log(clean_chunks))

  expect_identical(result$status, 0L)
  expect_contract_valid(result$path)
  expect_identical(result$document$result_kind, "package-check")
  expect_identical(result$document$status, "pass")
  expect_identical(
    result$document$check,
    list(
      errors = 0L,
      warnings = 0L,
      notes = 0L,
      log_path = "colocboost.Rcheck/00check.log"
    )
  )
})

test_that("the exact extension-type informational heading is accepted", {
  result <- run_parser(make_log(c(
    "* checking for file 'colocboost/DESCRIPTION' ... OK",
    "* checking extension type ... Package",
    "* checking tests ... INFO"
  )))

  expect_identical(result$status, 0L)
  expect_identical(result$document$status, "pass")
  expect_contract_valid(result$path)
})

test_that("a complete final Status line does not require a final newline", {
  result <- run_parser(make_log(clean_chunks), trailing_newline = FALSE)

  expect_identical(result$status, 0L)
  expect_identical(result$document$status, "pass")
})

test_that("ERROR and WARNING headings always fail with exact counts", {
  cases <- list(
    error = list(
      chunks = c("* checking examples ... ERROR", "actual error output"),
      footer = "1 ERROR",
      field = "errors"
    ),
    warning = list(
      chunks = c("* checking PDF version of manual ... WARNING", "warning output"),
      footer = "1 WARNING",
      field = "warnings"
    )
  )

  for (name in names(cases)) {
    fixture <- cases[[name]]
    result <- run_parser(make_log(fixture$chunks, fixture$footer))
    expect_identical(result$status, 1L, info = name)
    expect_contract_valid(result$path)
    expect_identical(result$document$result_kind, "package-check", info = name)
    expect_identical(result$document$status, "fail", info = name)
    expect_identical(result$document$check[[fixture$field]], 1L, info = name)
  }
})

test_that("unexpected NOTE fails and exact NOTE allowlist passes", {
  note_output <- make_note_rule()$output
  chunks <- c(
    "* checking installed package size ... NOTE",
    strsplit(note_output, "\n", fixed = TRUE)[[1L]]
  )

  unexpected <- run_parser(make_log(chunks, "1 NOTE"))
  exact <- run_parser(
    make_log(chunks, "1 NOTE"),
    policy = write_policy(list(make_note_rule()))
  )

  expect_identical(unexpected$status, 1L)
  expect_identical(unexpected$document$status, "fail")
  expect_identical(unexpected$document$check$notes, 1L)
  expect_contract_valid(unexpected$path)
  expect_identical(exact$status, 0L)
  expect_identical(exact$document$status, "pass")
  expect_identical(exact$document$check$notes, 1L)
  expect_contract_valid(exact$path)
})

test_that("NOTE matching is byte-exact apart from line endings", {
  rule <- make_note_rule()
  changed <- sub("5.0Mb", "5.1Mb", rule$output, fixed = TRUE)
  changed_whitespace <- sub("size is  ", "size is ", rule$output, fixed = TRUE)
  crlf <- gsub("\n", "\r\n", rule$output, fixed = TRUE)

  for (note_output in list(changed, changed_whitespace)) {
    result <- run_parser(
      make_log(c(
        "* checking installed package size ... NOTE",
        strsplit(note_output, "\n", fixed = TRUE)[[1L]]
      ), "1 NOTE"),
      policy = write_policy(list(rule))
    )
    expect_identical(result$status, 1L)
    expect_identical(result$document$status, "fail")
  }

  normalized <- run_parser(
    make_log(c(
      "* checking installed package size ... NOTE",
      strsplit(crlf, "\r\n", fixed = TRUE)[[1L]]
    ), "1 NOTE"),
    policy = write_policy(list(rule)),
    line_ending = "\r\n"
  )
  expect_identical(normalized$status, 0L)
})

test_that("unused, duplicate, expired, and multi-match NOTE rules fail closed", {
  rule <- make_note_rule()
  note_chunks <- c(
    "* checking installed package size ... NOTE",
    strsplit(rule$output, "\n", fixed = TRUE)[[1L]]
  )
  variants <- list(
    unused = run_parser(
      make_log(clean_chunks),
      policy = write_policy(list(rule))
    ),
    duplicate_observation = run_parser(
      make_log(c(note_chunks, note_chunks), "2 NOTEs"),
      policy = write_policy(list(rule))
    ),
    expired = run_parser(
      make_log(note_chunks, "1 NOTE"),
      policy = write_policy(list(make_note_rule(expires = "2000-01-01")))
    ),
    multi_match = run_parser(
      make_log(note_chunks, "1 NOTE"),
      policy = write_policy(list(rule, make_note_rule(id = "second-rule")))
    )
  )

  for (name in names(variants)) {
    result <- variants[[name]]
    expect_true(result$status != 0L, info = name)
    expect_contract_valid(result$path)
    expect_false(identical(result$document$status, "pass"), info = name)
  }
})

test_that("an ERROR cannot be converted into an allowed NOTE", {
  rule <- make_note_rule(
    check = "checking examples",
    output = "same output"
  )
  result <- run_parser(
    make_log(c("* checking examples ... ERROR", "same output"), "1 ERROR"),
    policy = write_policy(list(rule))
  )

  expect_identical(result$status, 1L)
  expect_identical(result$document$status, "fail")
  expect_identical(result$document$check$errors, 1L)
})

test_that("footer counts must exactly match stateful heading counts", {
  malformed <- list(
    count_mismatch = make_log(c("* checking examples ... NOTE", "x"), "2 NOTEs"),
    zero = make_log(clean_chunks, "0 NOTEs"),
    wrong_plural = make_log(c("* checking examples ... NOTE", "x"), "1 NOTEs"),
    wrong_order = make_log(
      c("* checking examples ... NOTE", "x", "* checking tests ... WARNING", "y"),
      "1 NOTE, 1 WARNING"
    ),
    duplicate_kind = make_log(
      c("* checking examples ... NOTE", "x", "* checking tests ... NOTE", "y"),
      "1 NOTE, 1 NOTE"
    ),
    unknown = make_log(clean_chunks, "1 FAILURE")
  )

  for (name in names(malformed)) {
    expect_infrastructure_failure(run_parser(malformed[[name]]))
  }
})

test_that("missing DONE, missing Status, missing headings, and truncated headings fail", {
  malformed <- list(
    missing_done = make_log(clean_chunks, done = FALSE),
    missing_status = make_log(clean_chunks, status = NULL),
    missing_heading = make_log(character()),
    truncated_heading = make_log("* checking tests ...")
  )

  for (name in names(malformed)) {
    expect_infrastructure_failure(run_parser(malformed[[name]]))
  }
})

test_that("duplicate footer markers and arbitrary footer garbage fail", {
  malformed <- list(
    duplicate_done = paste(
      make_log(clean_chunks, status = NULL),
      "* DONE",
      "Status: OK",
      sep = "\n"
    ),
    duplicate_status = paste(
      make_log(clean_chunks),
      "Status: OK",
      sep = "\n"
    ),
    trailing_garbage = paste(
      make_log(clean_chunks),
      "arbitrary trailing output",
      sep = "\n"
    )
  )

  for (name in names(malformed)) {
    expect_infrastructure_failure(run_parser(malformed[[name]]))
  }
})

test_that("Windows multi-arch parents and NONE or SKIPPED children are valid", {
  chunks <- c(
    "* checking examples ...",
    "** running examples for arch 'i386' ... OK",
    "** running examples for arch 'x64' ... NONE",
    "* checking tests ...",
    "** running tests for arch 'i386' ... SKIPPED",
    "** running tests for arch 'x64' ... OK",
    "* loading checks for arch 'x64'",
    "** checking whether package can be loaded ... OK"
  )
  result <- run_parser(make_log(chunks, elapsed = TRUE))

  expect_identical(result$status, 0L)
  expect_identical(result$document$status, "pass")
  expect_contract_valid(result$path)
})

test_that("malformed multi-arch parents fail and inner arch statuses are counted", {
  malformed_parent <- run_parser(make_log(c(
    "* checking examples ...",
    "parent output interrupts the required child adjacency",
    "** running examples for arch 'x64' ... OK"
  )))
  inner_warning <- run_parser(make_log(c(
    "* checking tests ...",
    "** running tests for arch 'i386' ... OK",
    "** running tests for arch 'x64' ... WARNING",
    "inner warning"
  ), "1 WARNING"))

  expect_infrastructure_failure(malformed_parent)
  expect_identical(inner_warning$status, 1L)
  expect_identical(inner_warning$document$check$warnings, 1L)
})

test_that("a nonzero check process exit fails even when the footer is OK", {
  result <- run_parser(make_log(clean_chunks), check_exit_code = 17L)

  expect_identical(result$status, 1L)
  expect_identical(result$document$result_kind, "package-check")
  expect_identical(result$document$status, "fail")
  expect_contract_valid(result$path)
})

test_that("malformed CLI writes infrastructure JSON when identity is recoverable", {
  result <- run_parser(make_log(clean_chunks), extra_args = "--unknown=value")
  expect_infrastructure_failure(result)
})
