#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- arguments[grepl("^--file=", arguments)]
if (length(file_argument) != 1L) stop("Cannot locate this contract test.")
test_file <- normalizePath(sub("^--file=", "", file_argument), mustWork = TRUE)
ci_directory <- normalizePath(file.path(dirname(test_file), ".."), mustWork = TRUE)
probe_path <- file.path(ci_directory, "probe-runtime.R")
source <- readLines(probe_path, warn = TRUE)
text <- paste(source, collapse = "\n")

required <- c(
  "lexical_r <-",
  "selected_r <- absolute(lexical_r",
  "selected_r_home <- system2(",
  'selected_r, "RHOME"',
  'active_r_home <- absolute(R.home()',
  "r_executable = lexical_r",
  "r_resolved = selected_r",
  "matrix_checksum <- sum(product)",
  'maps_path <- "/proc/self/maps"',
  'capabilities("long.double")',
  'config_value(selected_r, "CC")',
  'config_value(selected_r, "CXX")',
  'config_value(selected_r, "FC")',
  '"CHECK_ARGS"',
  '"VALGRIND_OPTS"'
)
missing <- required[!vapply(required, grepl, logical(1L), x = text, fixed = TRUE)]
if (length(missing)) {
  stop("Runtime probe contract is missing: ", paste(missing, collapse = ", "))
}
if (grepl("jsonlite", text, fixed = TRUE)) {
  stop("Runtime probe must not require an unprepared jsonlite installation.")
}
if (grepl('file.path(R.home("bin")', text, fixed = TRUE)) {
  stop("Runtime probe must compare R homes, not inequivalent launcher paths.")
}

function_lines <- source[
  seq.int(
    grep("^os_release_evidence <- function", source),
    grep("^write_evidence <- function", source) - 1L
  )
]
eval(parse(text = paste(function_lines, collapse = "\n")))

for (fixture in list(
  c("ID=ubuntu", 'VERSION_ID="24.04"'),
  c('ID="fedora"', "VERSION_ID=44")
)) {
  os_release <- tempfile("os-release-")
  on.exit(unlink(os_release), add = TRUE)
  writeLines(fixture, os_release, useBytes = TRUE)
  observed <- os_release_evidence(os_release)
  expected_id <- if (grepl("ubuntu", fixture[[1L]], fixed = TRUE)) "ubuntu" else "fedora"
  if (!identical(observed$id, expected_id)) stop("OS release ID was not parsed exactly.")
}

cat("Task 10 runtime policy contract passed.\n")
