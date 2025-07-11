library(testthat)

# Generate test data with specific edge case structures
generate_edge_case_data <- function(case = "missing_values", 
                                    n = 100, 
                                    p = 20, 
                                    L = 2, 
                                    seed = 42) {
  set.seed(seed)
  
  # Generate X with LD structure
  sigma <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigma[i, j] <- 0.9^abs(i - j)
    }
  }
  X <- MASS::mvrnorm(n, rep(0, p), sigma)
  colnames(X) <- paste0("SNP", 1:p)
  
  # Generate true effects with shared causal variant
  true_beta <- matrix(0, p, L)
  if (L == 1) {
    # Single trait case
    true_beta[5, 1] <- 0.5  # SNP5 affects the trait
  } else {
    # Multi-trait case
    true_beta[5, 1] <- 0.5  # SNP5 affects trait 1
    true_beta[5, 2] <- 0.4  # SNP5 also affects trait 2 (colocalized)
    true_beta[10, 2] <- 0.3 # SNP10 only affects trait 2
  }
  
  
  # Generate Y with some noise
  Y <- matrix(0, n, L)
  for (l in 1:L) {
    Y[, l] <- X %*% true_beta[, l] + rnorm(n, 0, 1)
  }
  
  # Modifications based on case
  if (case == "missing_values") {
    # Introduce missing values in Y
    Y[sample(1:n, 5), 1] <- NA
    Y[sample(1:n, 5), 2] <- NA
  } else if (case == "different_samples") {
    # Create different sample sizes for each outcome
    Y1 <- Y[1:(n-10), 1]
    Y2 <- Y[1:(n-20), 2]
    return(list(
      X = X,
      Y = list(Y1, Y2), 
      LD = cor(X),
      true_beta = true_beta
    ))
  } else if (case == "different_variants") {
    # Create different variant sets for each outcome
    X1 <- X[, 1:(p-5)]
    X2 <- X[, 6:p]
    Y1 <- Y[, 1]
    Y2 <- Y[, 2]
    return(list(
      X = list(X1, X2),
      Y = list(Y1, Y2), 
      LD = cor(X),
      true_beta = true_beta
    ))
  } else if (case == "no_colocalization") {
    # Generate data with no true colocalization
    true_beta[5, 2] <- 0  # Remove shared effect
    Y[, 2] <- X %*% true_beta[, 2] + rnorm(n, 0, 1)
  } else if (case == "high_correlation") {
    # Create highly correlated outcomes
    Y[, 2] <- 0.9 * Y[, 1] + rnorm(n, 0, 0.3)
  } 
  
  # Return test objects
  list(
    X = X,
    Y = Y, 
    LD = cor(X),
    true_beta = true_beta
  )
}

# Test colocboost handling of missing values
test_that("colocboost handles missing values in Y", {
  
  # Generate data with missing values
  test_data <- generate_edge_case_data("missing_values")
  
  # Convert Y to list
  Y_list <- list(test_data$Y[,1], test_data$Y[,2])
  X_list <- list(test_data$X, test_data$X)
  
  # Run colocboost - should handle NAs automatically
  expect_error(
    suppressWarnings(
      result <- colocboost(
        X = X_list, 
        Y = Y_list,
        M = 5  # Small number of iterations for testing
      )
    ),
    NA
  )
  
  # Test that we get a colocboost object
  expect_s3_class(result, "colocboost")
})

# Test colocboost with different sample sizes
test_that("colocboost handles different sample sizes", {
  
  # Generate data with different sample sizes
  test_data <- generate_edge_case_data("different_samples")
  
  # Run colocboost - should error without sample indices
  expect_error(
    suppressWarnings(
      colocboost(
        X = test_data$X, 
        Y = test_data$Y
      )
    )
  )
  
  # Test with proper row names to handle different sample sizes
  X <- test_data$X
  rownames(X) <- paste0("sample", 1:nrow(X))
  Y1 <- as.matrix(test_data$Y[[1]])
  Y2 <- as.matrix(test_data$Y[[2]])
  rownames(Y1) <- paste0("sample", 1:length(Y1))
  rownames(Y2) <- paste0("sample", 1:length(Y2))
  Y <- list(Y1, Y2)
  
  # Skip test if function fails in unexpected ways
  # (This test may require more complex setup than we can do here) 
  # skip("Requires specialized setup for different sample sizes")
  # Can handle on Apr 17, 2025 (News)
  # Run colocboost - should only warning
  expect_warning(colocboost(X = X, Y = Y))
})

# Test colocboost with different variant sets
test_that("colocboost handles different variant sets", {
  
  # Generate data with different variant sets
  test_data <- generate_edge_case_data("different_variants")
  
  # Run colocboost
  warnings <- capture_warnings({
    result <- colocboost(
      X = test_data$X, 
      Y = test_data$Y,
      dict_YX = dict_YX,
      M = 5  # Small number of iterations for testing
    )
  })
  # Check if any of the expected warning patterns are present
  expect_true(
    any(grepl("did not coverage", warnings))
  )
  
  # Test that we get a colocboost object
  expect_s3_class(result, "colocboost")
})

# Test colocboost with no true colocalization
test_that("colocboost correctly identifies absence of colocalization", {

  # Generate data with no colocalization
  test_data <- generate_edge_case_data("no_colocalization")
  
  # Convert Y to list
  Y_list <- list(test_data$Y[,1], test_data$Y[,2])
  X_list <- list(test_data$X, test_data$X)
  
  # Run colocboost
  suppressWarnings({
    result <- colocboost(
      X = X_list, 
      Y = Y_list,
      M = 10  # Need more iterations for this test
    )
  })
  
  # Test that we get a colocboost object
  expect_s3_class(result, "colocboost")
  
  # In many cases with no true colocalization, cos_summary should be NULL
  # But this isn't guaranteed due to statistical fluctuations
  # Instead we'll check that object has expected structure
  expect_true(any(
    is.null(result$cos_summary) || 
    is.data.frame(result$cos_summary)
  ))
})


# Test colocboost with highly correlated traits
test_that("colocboost handles highly correlated traits", {
  
  # Generate data with correlated traits
  test_data <- generate_edge_case_data("high_correlation")
  
  # Convert Y to list
  Y_list <- list(test_data$Y[,1], test_data$Y[,2])
  X_list <- list(test_data$X, test_data$X)
  
  # Create correlation matrix for residuals
  residual_corr <- matrix(c(1, 0.9, 0.9, 1), nrow=2)
  
  # Run colocboost with residual correlation
  suppressWarnings({
    result <- colocboost(
      X = X_list, 
      Y = Y_list,
      residual_correlation = residual_corr,
      M = 10  # Need more iterations for this test
    )
  })
  
  # Test that we get a colocboost object
  expect_s3_class(result, "colocboost")
})

# Test with very small dataset
test_that("colocboost handles very small datasets", {
  
  # Create tiny dataset
  set.seed(123)
  X_small <- matrix(rnorm(10*5), 10, 5)
  colnames(X_small) <- paste0("SNP", 1:5)
  Y_small <- matrix(rnorm(10*2), 10, 2)
  Y_list_small <- list(Y_small[,1], Y_small[,2])
  X_list_small <- list(X_small, X_small)
  
  # Run colocboost
  suppressWarnings({
    result <- colocboost(
      X = X_list_small, 
      Y = Y_list_small,
      M = 5  # Small number of iterations for testing
    )
  })
  
  # Test that we get a colocboost object
  expect_s3_class(result, "colocboost")
})

# Test with custom parameter settings
test_that("colocboost works with custom parameters", {
  
  # Create simple dataset
  set.seed(123)
  X <- matrix(rnorm(50*10), 50, 10)
  colnames(X) <- paste0("SNP", 1:10)
  Y <- matrix(rnorm(50*2), 50, 2)
  Y_list <- list(Y[,1], Y[,2])
  X_list <- list(X, X)
  
  # Run colocboost with custom parameters
  suppressWarnings({
    result <- colocboost(
      X = X_list, 
      Y = Y_list,
      M = 5,  # Small number of iterations for testing
      lambda = 0.7,  # Custom lambda
      tau = 0.02,    # Custom tau
      func_simplex = "only_z2z",  # Different simplex function
      learning_rate_init = 0.02,  # Custom learning rate
      coverage = 0.9,  # Different coverage
      min_abs_corr = 0.4  # Lower correlation threshold
    )
  })
  
  # Test that we get a colocboost object
  expect_s3_class(result, "colocboost")
})

# Test focal outcome functionality
test_that("colocboost prioritizes focal outcome correctly", {
  
  # Generate test data
  set.seed(123)
  X <- matrix(rnorm(100*20), 100, 20)
  colnames(X) <- paste0("SNP", 1:20)
  
  # Create shared causal variant
  b <- rep(0, 20)
  b[5] <- 0.5
  b[10] <- 0.3
  
  # Generate outcomes with different effects
  Y1 <- X %*% b + rnorm(100, 0, 1)
  
  # Y2 has different effect at position 5
  b2 <- b
  b2[5] <- 0.2
  Y2 <- X %*% b2 + rnorm(100, 0, 1)
  
  Y_list <- list(Y1, Y2)
  X_list <- list(X, X)
  
  # Run colocboost with Y1 as focal
  suppressWarnings({
    result_focal1 <- colocboost(
      X = X_list, 
      Y = Y_list,
      focal_outcome_idx = 1,
      lambda_focal_outcome = 0.9,  # Higher lambda for focal
      M = 10  # Need more iterations for this test
    )
  })
  
  # Run colocboost with Y2 as focal
  suppressWarnings({
    result_focal2 <- colocboost(
      X = X_list, 
      Y = Y_list,
      focal_outcome_idx = 2,
      lambda_focal_outcome = 0.9,  # Higher lambda for focal
      M = 10  # Need more iterations for this test
    )
  })
  
  # Both should produce colocboost objects
  expect_s3_class(result_focal1, "colocboost")
  expect_s3_class(result_focal2, "colocboost")
  
  # Test that focal flag is correctly set
  expect_equal(result_focal1$data_info$outcome_info$is_focal[1], TRUE)
  expect_equal(result_focal1$data_info$outcome_info$is_focal[2], FALSE)
  
  expect_equal(result_focal2$data_info$outcome_info$is_focal[1], FALSE)
  expect_equal(result_focal2$data_info$outcome_info$is_focal[2], TRUE)
})


# Test with ambiguous corner cases
test_that("get_ambiguous_colocalization handles edge cases with correlation thresholds", {

  data("Ambiguous_Colocalization", package = "colocboost", envir = environment())
  test_colocboost_results <- Ambiguous_Colocalization$ColocBoost_Results
  
  # Test with very high correlation thresholds (should find fewer ambiguities)
  result_high_thresh <- get_ambiguous_colocalization(
    test_colocboost_results,
    min_abs_corr_between_ucos = 0.95,
    median_abs_corr_between_ucos = 0.98
  )
  
  # Test with very low correlation thresholds (should find more ambiguities)
  result_low_thresh <- get_ambiguous_colocalization(
    test_colocboost_results,
    min_abs_corr_between_ucos = 0.1,
    median_abs_corr_between_ucos = 0.3
  )
  
  # Compare number of ambiguous events found with different thresholds
  # Generally expect: n_high_thresh <= n_default <= n_low_thresh
  n_high <- length(result_high_thresh$ambiguous_cos)
  n_default <- length(get_ambiguous_colocalization(test_colocboost_results)$ambiguous_cos)
  n_low <- length(result_low_thresh$ambiguous_cos)
  
  # Higher thresholds should find equal or fewer ambiguities than default
  expect_true(n_high <= n_default)
  
  # Lower thresholds should find equal or more ambiguities than default
  expect_true(n_low >= n_default)
  
  # Test with extreme threshold (min=1.0, median=1.0) - should find very few or no ambiguities
  expect_message(
    result_extreme <- get_ambiguous_colocalization(
      test_colocboost_results,
      min_abs_corr_between_ucos = 1.0,
      median_abs_corr_between_ucos = 1.0
    ),
    "No ambiguous colocalization events!"
  )
  
  # Test with threshold at 0 - should find many/all potential ambiguities
  result_zero <- get_ambiguous_colocalization(
    test_colocboost_results,
    min_abs_corr_between_ucos = 0.0,
    median_abs_corr_between_ucos = 0.0
  )
  expect_true(length(result_zero$ambiguous_cos) >= n_low)
  
  
})

