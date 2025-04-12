setwd("../")

# Install necessary packages
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools not found")
}
if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("testthat not found")
}
if (!requireNamespace("covr", quietly = TRUE)) {
  stop("covr not found")
}
if (!requireNamespace("roxygen2", quietly = TRUE)) {
  stop("roxygen2 not found")
}

# Set up testthat infrastructure
if (!dir.exists("tests/testthat")) {
  devtools::use_testthat()
}

# Create GitHub Actions workflow for continuous integration
if (!dir.exists(".github/workflows")) {
  dir.create(".github/workflows", recursive = TRUE, showWarnings = FALSE)
}

# Write GitHub Actions workflow file
workflow_file <- ".github/workflows/R-CMD-check.yaml"
writeLines(
'
on:
  push:
    branches: [main, master]
  pull_request:
    branches: [main, master]

name: R-CMD-check

jobs:
  R-CMD-check:
    runs-on: ${{ matrix.config.os }}

    name: ${{ matrix.config.os }} (${{ matrix.config.r }})

    strategy:
      fail-fast: false
      matrix:
        config:
          - {os: windows-latest, r: "release"}
          - {os: macOS-latest, r: "release"}
          - {os: ubuntu-latest, r: "release"}

    env:
      R_REMOTES_NO_ERRORS_FROM_WARNINGS: true
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}

    steps:
      - uses: actions/checkout@v2

      - uses: r-lib/actions/setup-r@v2
        with:
          r-version: ${{ matrix.config.r }}

      - uses: r-lib/actions/setup-pandoc@v2

      - name: Install dependencies
        run: |
          install.packages(c("remotes", "rcmdcheck", "covr", "testthat"))
          remotes::install_deps(dependencies = TRUE)
        shell: Rscript {0}

      - name: Check
        run: rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"), error_on = "warning", check_dir = "check")
        shell: Rscript {0}

      - name: Test coverage
        run: covr::codecov()
        shell: Rscript {0}
', workflow_file)

# Create a test coverage workflow
coverage_file <- ".github/workflows/test-coverage.yaml"
writeLines(
'
on:
  push:
    branches: [main, master]
  pull_request:
    branches: [main, master]

name: test-coverage

jobs:
  test-coverage:
    runs-on: ubuntu-latest
    
    steps:
      - uses: actions/checkout@v2
      
      - uses: r-lib/actions/setup-r@v2
      
      - name: Install dependencies
        run: |
          install.packages(c("remotes", "covr"))
          remotes::install_deps(dependencies = TRUE)
        shell: Rscript {0}
        
      - name: Test coverage
        run: covr::codecov()
        shell: Rscript {0}
', coverage_file)

# Add code coverage badge to README.md
readme_file <- "README.md"
if (file.exists(readme_file)) {
  readme_content <- readLines(readme_file)
  if (!any(grepl("codecov", readme_content))) {
    badge <- "[![Codecov test coverage](https://codecov.io/gh/StatFunGen/colocboost/branch/master/graph/badge.svg)](https://codecov.io/gh/StatFunGen/colocboost?branch=master)"
    # Add badge after the first line if it's a title
    if (length(readme_content) > 0) {
      readme_content <- c(readme_content[1], badge, readme_content[-1])
    } else {
      readme_content <- c("# colocboost", badge)
    }
    writeLines(readme_content, readme_file)
  }
}

# Set up basic testthat structure
message("testthat setup complete. Next, create test files for each R file in the package.")