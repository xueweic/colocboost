library(testthat)

# Utility function to generate test data
generate_test_data <- function(n = 100, p = 20, L = 2, seed = 42) {
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
  
  # Generate true effects - create a shared causal variant
  true_beta <- matrix(0, p, L)
  true_beta[5, 1] <- 0.5  # SNP5 affects trait 1
  true_beta[5, 2] <- 0.4  # SNP5 also affects trait 2 (colocalized)
  true_beta[10, 2] <- 0.3 # SNP10 only affects trait 2
  
  # Generate Y with some noise
  Y <- matrix(0, n, L)
  for (l in 1:L) {
    Y[, l] <- X %*% true_beta[, l] + rnorm(n, 0, 1)
  }
  
  # Prepare LD matrix
  LD <- cor(X)
  
  # Return test objects
  list(
    X = X,
    Y = Y, 
    LD = LD,
    true_beta = true_beta
  )
}

# Create small dataset for basic tests
test_data <- generate_test_data()

# Basic test for colocboost functionality with individual data
test_that("colocboost runs with individual data", {
  
  # Convert Y to list
  Y_list <- list(test_data$Y[,1], test_data$Y[,2])
  X_list <- list(test_data$X, test_data$X)
  
  # Run colocboost with minimal parameters
  warnings <- capture_warnings({
    result <- colocboost(
      X = X_list, 
      Y = Y_list,
      M = 10,  # Small number of iterations for testing
      output_level = 2  # More detailed output for testing
    )
  })

  # Check if any of the expected warning patterns are present
  expect_true(
    any(grepl("smallest number of variables", warnings)) || 
    any(grepl("did not coverage", warnings))
  )
  
  # Test that we get a colocboost object
  expect_s3_class(result, "colocboost")
  
  # Check structure of results
  expect_type(result$data_info, "list")
  expect_type(result$model_info, "list")
  
  # Check dimensions
  expect_equal(length(result$data_info$variables), ncol(test_data$X))
  expect_equal(result$data_info$n_outcomes, 2)
})

# Test with summary statistics
test_that("colocboost runs with summary statistics", {
  
  # Generate summary statistics from the individual data
  X <- test_data$X
  Y <- test_data$Y
  LD <- test_data$LD
  
  # Calculate beta, se, and z-scores
  beta <- matrix(0, ncol(X), ncol(Y))
  se <- matrix(0, ncol(X), ncol(Y))
  z <- matrix(0, ncol(X), ncol(Y))
  
  for (i in 1:ncol(Y)) {
    fit <- lm(Y[,i] ~ X)
    beta[,i] <- coef(fit)[-1]
    se[,i] <- summary(fit)$coefficients[-1, "Std. Error"]
    z[,i] <- beta[,i]/se[,i]
  }
  
  # Create summary statistics data frames
  sumstat_list <- list()
  for (i in 1:ncol(Y)) {
    sumstat_list[[i]] <- data.frame(
      beta = beta[,i],
      sebeta = se[,i],
      z = z[,i],
      n = nrow(X),
      variant = colnames(X)
    )
  }
  
  # Run colocboost with summary statistics
  warnings <- capture_warnings({
    result <- colocboost(
      sumstat = sumstat_list,
      LD = list(LD),
      M = 10,  # Small number of iterations for testing
      output_level = 2  # More detailed output for testing
    )
  })

  # Check if any of the expected warning patterns are present
  expect_true(
    any(grepl("smallest number of variables", warnings)) || 
    any(grepl("did not coverage", warnings))
  )
  
  # Test that we get a colocboost object
  expect_s3_class(result, "colocboost")
  
  # Check dimensions
  expect_equal(length(result$data_info$variables), ncol(X))
  expect_equal(result$data_info$n_outcomes, 2)
})



# Test get_cos_summary functionality
test_that("get_cos_summary returns expected structure", {
  
  # Convert Y to list
  Y_list <- list(test_data$Y[,1], test_data$Y[,2])
  X_list <- list(test_data$X, test_data$X)
  
  # Run colocboost with minimal parameters
  warnings <- capture_warnings({
    result <- colocboost(
      X = X_list, 
      Y = Y_list,
      M = 10,  # Small number of iterations for testing
      output_level = 2  # More detailed output for testing
    )
  })

  # Check if any of the expected warning patterns are present
  expect_true(
    any(grepl("smallest number of variables", warnings)) || 
    any(grepl("did not coverage", warnings))
  )
  
  # Get summary
  if (!is.null(result$cos_summary)) {
    summary <- result$cos_summary
    
    # Check structure of summary table
    expect_true(is.data.frame(summary))
    
    # Check expected columns exist
    expected_cols <- c("focal_outcome", "colocalized_outcomes", "cos_id", 
                      "purity", "top_variable", "top_variable_vcp")
    for (col in expected_cols) {
      expect_true(col %in% colnames(summary))
    }
  } else {
    # If no colocalization is found, just check that the summary is NULL
    expect_null(result$cos_summary)
  }
})


# Test focal outcome functionality
test_that("colocboost handles focal outcome correctly", {
  
  # Convert Y to list
  Y_list <- list(test_data$Y[,1], test_data$Y[,2])
  X_list <- list(test_data$X, test_data$X)
  
  # Run colocboost with focal_outcome_idx = 1
  warnings <- capture_warnings({
    result <- colocboost(
      X = X_list, 
      Y = Y_list,
      focal_outcome_idx = 1,
      M = 10,  # Small number of iterations for testing
      output_level = 2  # More detailed output for testing
    )
  })

  # Check if any of the expected warning patterns are present
  expect_true(
    any(grepl("smallest number of variables", warnings)) || 
    any(grepl("did not coverage", warnings))
  )
  
  # Test that we get a colocboost object
  expect_s3_class(result, "colocboost")
  
  # Check focal outcome is correctly set
  expect_equal(result$data_info$outcome_info$is_focal[1], TRUE)
  expect_equal(result$data_info$outcome_info$is_focal[2], FALSE)
})

# Test colocboost_plot functionality (basic call should not error)
test_that("colocboost_plot runs without error", {
  
  # Convert Y to list
  Y_list <- list(test_data$Y[,1], test_data$Y[,2])
  X_list <- list(test_data$X, test_data$X)
  
  # Run colocboost with minimal parameters
  warnings <- capture_warnings({
    result <- colocboost(
      X = X_list, 
      Y = Y_list,
      M = 10,  # Small number of iterations for testing
      output_level = 2  # More detailed output for testing
    )
  })

  # Check if any of the expected warning patterns are present
  expect_true(
    any(grepl("smallest number of variables", warnings)) || 
    any(grepl("did not coverage", warnings))
  )
  
  # Basic plotting should not throw an error
  # Use suppressWarnings as the function might warn about lack of colocalization
  expect_error(suppressWarnings(colocboost_plot(result)), NA)
})

# Test get_robust_colocalization functionality
test_that("get_robust_colocalization maintains colocboost structure", {
  
  # Convert Y to list
  Y_list <- list(test_data$Y[,1], test_data$Y[,2])
  X_list <- list(test_data$X, test_data$X)
  
  # Run colocboost with minimal parameters
  warnings <- capture_warnings({
    result <- colocboost(
      X = X_list, 
      Y = Y_list,
      M = 10,
      output_level = 2
    )
  })

  # Check if any of the expected warning patterns are present
  expect_true(
    any(grepl("smallest number of variables", warnings)) || 
    any(grepl("did not coverage", warnings))
  )
  
  # Run get_strong_colocalization
  strong_result <- get_robust_colocalization(
    result, 
    cos_npc_cutoff = 0.2,  # Lower threshold for testing
    npc_outcome_cutoff = 0.1  # Lower threshold for testing
  )
  
  # Should still be a colocboost object
  expect_s3_class(strong_result, "colocboost")
  
  # Should have the same data_info
  expect_equal(strong_result$data_info, result$data_info)
})

# Test error handling with malformed input
test_that("colocboost handles missing/invalid inputs appropriately", {

  # Test missing both individual data and summary stats
  expect_warning(colocboost(), "Error: No individual data")
  
  # Test mismatched dimensions
  X_bad <- test_data$X[1:(nrow(test_data$X) - 10), ]
  Y_list <- list(test_data$Y[,1], test_data$Y[,2])
  X_list <- list(X_bad, X_bad)
  
  expect_error(
    suppressWarnings(
      colocboost(
        X = X_list, 
        Y = Y_list
      )
    ),
    "Please provide the sample index of X and Y, since they do not have the same samples!"
  )
})


