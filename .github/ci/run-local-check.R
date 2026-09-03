#!/usr/bin/env Rscript

devtools::check(
  pkg = ".",
  args = "--as-cran",
  build_args = character(),
  manual = TRUE
)
