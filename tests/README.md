# colocboost unit-testing

## Overview

This repository contains a comprehensive testing framework for the [colocboost](https://github.com/StatFunGen/colocboost) R package. The framework is designed to ensure the reliability and correctness of package functionality through automated testing.

## Quick Start
   
1. Our unit testing setup is managed by `pixi`.  Please follow the instructions at [https://pixi.sh/latest/#installation](https://pixi.sh/latest/#installation) if you have not already installed `pixi`.

2. In the root of this repository, run the helper script to create a pixi.toml file.  This file is deliberately ignored in `.gitignore` because it is ephemeral and will be regenerated whenever CI is run.

    ```bash
    .github/workflows/create_toml_from_yaml.sh $(pwd)
    ```

3. Run all tests with a `pixi` task:
   ```bash
   pixi run devtools_test
   ```
   or test one file using `pixi run`:
   ```bash
   pixi run R -e 'devtools::test_active_file("tests/testthat/test_colocboost.R")'
   ```

## Files and Structure

- `testthat/`: Directory containing test files
  - `test_package.R`: Tests for basic package functionality
  - `test_colocboost.R`: Tests for the main colocalization functions
  - `test_inference.R`: Tests for the post inference functions
  - `test_utils.R`: Tests for utility functions
  - `test_model.R`: Tests for model fitting and prediction functions
- `.github/workflows/`: GitHub Actions workflow configurations

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

The `covr` package is used to measure test coverage, which indicates what percentage of the code is being tested.

To generate a test coverage report:
```bash
pixi run codecov
```

## GitHub Actions Workflow

This testing framework includes GitHub Actions workflows that automatically run the tests on every pull request. The workflows test the package on Linux and macOS to ensure cross-platform compatibility.

The workflow is defined in `.github/workflows/ci.yaml` and automatically runs:
- `R CMD check`
- Test coverage reporting with [codecov.io](codecov.io)

Test results and coverage statistics are available on the GitHub Actions page after each pull request.

For more information, check the [testthat documentation](https://testthat.r-lib.org/).
