#!/usr/bin/env Rscript

# Script to install the package and run tests

# Check for required packages
required_packages <- c("devtools", "testthat", "covr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Cannot find required package: %s\n", pkg))
  }
}

# Path to working directory
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  work_dir <- args[1]
} else {
  work_dir <- getwd()
}

setwd(work_dir)
cat(sprintf("Working directory: %s\n", work_dir))

# Install the package
# cat("Installing colocboost package...\n")
# devtools::install(".", dependencies = FALSE, quiet = TRUE)

# Run tests
cat("Running tests...\n")
devtools::load_all('../')
testthat::test_dir("testthat/")

# Calculate test coverage
cat("Calculating test coverage...\n")
coverage <- covr::package_coverage()
print(coverage)
covr::report(coverage)

cat("Tests completed.\n")