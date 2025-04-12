# colocboost unit-testing

## Overview

This repository contains a comprehensive testing framework for the [colocboost](https://github.com/StatFunGen/colocboost) R package. The framework is designed to ensure the reliability and correctness of the package's functionality through automated testing.

## Quick Start
   
1. Navigate to test folder:
   ```bash
   cd tests
   ```

2. **First time use**: run the setup script to install required packages and configure the testing environment:
   ```r
   # install.packages(c("devtools", "testthat", "covr", "roxygen2"))
   source("setup_testthat.R")
   ```

3. Run all tests:
   ```r
   devtools::load_all()
   devtools::test()
   ```  
   or,
   ```bash
   Rscript run_tests.R
   ```
   To test one file:
   ```r
   devtools::test_active_file("testthat/test_colocboost.R")
   ```

## Files and Structure

- `setup_testthat.R`: Script to set up the testthat infrastructure
- `run_tests.R`: Script to run all tests and generate test coverage reports
- `testthat/`: Directory containing test files
  - `test_package.R`: Tests for basic package functionality
  - `test_colocboost.R`: Tests for the main colocalization functions
  - `test_utils.R`: Tests for utility functions
  - `test_model.R`: Tests for model fitting and prediction functions
- `.github/workflows/`: GitHub Actions workflow configurations

## How To Use

### Running Tests Locally

To run the tests locally, you can use:

```r
devtools::test()
```

Or run individual test files:

```r
devtools::test_file("testthat/test_colocboost.R")
```

### Adding New Tests

When adding new functionality to the colocboost package, corresponding tests should be added to maintain test coverage. Follow these steps:

1. Identify which test file should contain the new tests, or create a new test file if necessary
2. Write test functions using the testthat package's expectations
3. Run the tests to ensure they pass

Example test:

```r
test_that("new_function produces expected output", {
  # Arrange
  input_data <- prepare_test_data()
  
  # Act
  result <- new_function(input_data)
  
  # Assert
  expect_equal(result$some_value, expected_value)
  expect_true(is.data.frame(result$data))
})
```

### Test Coverage

The `covr` package is used to measure test coverage, which indicates what percentage of your code is being tested. Aim for at least 80% coverage for a reliable package.

To generate a test coverage report:
```r
library(covr)
coverage <- package_coverage()
report(coverage)
```

## GitHub Actions Workflow

This testing framework includes GitHub Actions workflows that automatically run the tests on every push and pull request. The workflows test the package on multiple operating systems (Windows, macOS, and Linux) to ensure cross-platform compatibility.

The workflow is defined in `.github/workflows/R-CMD-check.yaml` and automatically runs:
- R CMD check
- Test coverage reporting

Test results and coverage statistics are available on the GitHub Actions page after each push or pull request.

For more information, check the [testthat documentation](https://testthat.r-lib.org/).
