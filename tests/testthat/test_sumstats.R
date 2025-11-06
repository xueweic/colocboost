library(testthat)

# Copy generate_test_data from test_colocboost.R
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

# Use the existing test data generation function
# But extend it to create summary statistics directly
generate_sumstat_test_data <- function(n = 100, p = 20, L = 2, seed = 42) {
  # Get individual level data first
  test_data <- generate_test_data(n, p, L, seed)
  
  # Generate summary statistics from individual data
  X <- test_data$X
  Y <- test_data$Y
  LD <- test_data$LD
  
  # Calculate beta, se, and z-scores
  beta <- matrix(0, ncol(X), ncol(Y))
  se <- matrix(0, ncol(X), ncol(Y))
  z <- matrix(0, ncol(X), ncol(Y))
  
  for (i in 1:ncol(Y)) {
    output = matrix(0,ncol(X),2)
    for (j in 1:ncol(X)) {
      fit = summary(lm(Y[,i] ~ X[,j]))$coef
      if (nrow(fit) == 2)
        output[j,] = as.vector(fit[2,1:2])
      else
        output[j,] = c(0,0)
    }
    beta[,i] <- output[,1]
    se[,i] <- output[,2]
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
  
  # Add effect matrices for HyPrColoc format
  rownames(beta) <- rownames(se) <- colnames(X)
  colnames(beta) <- colnames(se) <- paste0("Y", 1:ncol(Y))
  
  # Return all formats
  list(
    X = X,
    Y = Y,
    LD = LD,
    sumstat = sumstat_list,
    effect_est = beta,
    effect_se = se,
    effect_n = rep(n, ncol(Y))
  )
}

# Create summary statistics test data
test_sumstat_data <- generate_sumstat_test_data()

# Test 1: Basic summary statistics input
test_that("colocboost runs with basic summary statistics format", {
  # Run colocboost with sumstat and single LD matrix
  warnings <- capture_warnings({
    result <- colocboost(
        sumstat = test_sumstat_data$sumstat,
        LD = test_sumstat_data$LD,
        M = 10,  # Small number of iterations for testing
        output_level = 2  # More detailed output for testing
    )
  })

  # Check if any of the expected warning patterns are present
  expect_true(
    any(grepl("did not coverage", warnings))
  )
  
  # Test that we get a colocboost object
  expect_s3_class(result, "colocboost")
  
  # Check structure of results
  expect_type(result$data_info, "list")
  expect_type(result$model_info, "list")
  
  # Check dimensions
  expect_equal(length(result$data_info$variables), ncol(test_sumstat_data$X))
  expect_equal(result$data_info$n_outcomes, 2)
})

# Test 2: HyPrColoc compatible format
test_that("colocboost runs with HyPrColoc compatible format", {
  # Run colocboost with effect matrices
  warnings <- capture_warnings({
    result <- colocboost(
        effect_est = test_sumstat_data$effect_est,
        effect_se = test_sumstat_data$effect_se,
        effect_n = test_sumstat_data$effect_n,
        LD = test_sumstat_data$LD,
        M = 10,  # Small number of iterations for testing
        output_level = 2  # More detailed output for testing
    )
  })

  # Check if any of the expected warning patterns are present
  expect_true(
    any(grepl("did not coverage", warnings))
  )
  
  # Test that we get a colocboost object
  expect_s3_class(result, "colocboost")
  
  # Check dimensions
  expect_equal(length(result$data_info$variables), ncol(test_sumstat_data$X))
  expect_equal(result$data_info$n_outcomes, 2)
})

# Test 3: Multiple LD matrices with summary statistics
test_that("colocboost runs with matched LD matrices", {
  # Create a list of identical LD matrices for demonstration
  LD_list <- list(test_sumstat_data$LD, test_sumstat_data$LD)
  
  # Run colocboost with summary statistics and multiple LD matrices
  warnings <- capture_warnings({
    result <- colocboost(
        sumstat = test_sumstat_data$sumstat,
        LD = LD_list,
        M = 10,  # Small number of iterations for testing
        output_level = 2  # More detailed output for testing
    )
  })

  # Check if any of the expected warning patterns are present
  expect_true(
    any(grepl("did not coverage", warnings))
  )
  
  # Test that we get a colocboost object
  expect_s3_class(result, "colocboost")
  
  # Check dimensions
  expect_equal(length(result$data_info$variables), ncol(test_sumstat_data$X))
  expect_equal(result$data_info$n_outcomes, 2)
})

# Test 4: Arbitrary LD matrices with dictionary mapping
test_that("colocboost runs with dictionary-mapped LD matrices", {
  # Create a list of identical LD matrices for demonstration
  LD_list <- list(test_sumstat_data$LD, test_sumstat_data$LD)
  
  # Create a dictionary mapping sumstat indices to LD indices
  # First sumstat uses first LD, second sumstat uses second LD
  dict_sumstatLD <- matrix(c(1:2, 1:2), ncol = 2)
  
  # Run colocboost with dictionary mapping
  warnings <- capture_warnings({
    result <- colocboost(
        sumstat = test_sumstat_data$sumstat,
        LD = LD_list,
        dict_sumstatLD = dict_sumstatLD,
        M = 10,  # Small number of iterations for testing
        output_level = 2  # More detailed output for testing
    )
  })

  # Check if any of the expected warning patterns are present
  expect_true(
    any(grepl("did not coverage", warnings))
  )
  
  # Test that we get a colocboost object
  expect_s3_class(result, "colocboost")
  
  # Check dimensions
  expect_equal(length(result$data_info$variables), ncol(test_sumstat_data$X))
  expect_equal(result$data_info$n_outcomes, 2)
})

# Test 5: Missing LD matrix (LD-free approach)
test_that("colocboost runs with missing LD matrix", {
  # Run colocboost without LD matrix (should default to diagonal)
  warnings <- capture_warnings({
    result <- colocboost(
        sumstat = test_sumstat_data$sumstat,
        M = 1,  # Only one iteration expected in LD-free case
        output_level = 2  # More detailed output for testing
    )
  })

  # Check if any of the expected warning patterns are present
  expect_true(
    any(grepl("The smallest number of variables across outcomes is 20 < 100", warnings)) 
  )
  
  # Test that we get a colocboost object
  expect_s3_class(result, "colocboost")
  
  # Check dimensions
  expect_equal(length(result$data_info$variables), ncol(test_sumstat_data$X))
  expect_equal(result$data_info$n_outcomes, 2)
})

# Test 6: Error handling for mismatched inputs
test_that("colocboost handles mismatched inputs correctly", {
  # Create mismatched effect matrices
  bad_effect_est <- test_sumstat_data$effect_est[1:(nrow(test_sumstat_data$effect_est)-1), ]
  
  # Expect error with mismatched dimensions
  expect_warning(
    colocboost_validate_input_data(
      effect_est = bad_effect_est,
      effect_se = test_sumstat_data$effect_se
    ),
    "Error: effect_est and effect_se should be the same dimension!"
  )
})

# Test 7: Summary statistics without sample size
test_that("colocboost handles missing sample size correctly", {
  # Create summary stats without sample size
  sumstat_no_n <- lapply(test_sumstat_data$sumstat, function(df) {
    df$n <- NULL
    return(df)
  })
  
  # Should run but give a warning about missing sample size
  warnings <- capture_warnings({
    result <- colocboost(
      sumstat = sumstat_no_n,
      LD = test_sumstat_data$LD,
      M = 10,
      output_level = 2
    )
  })

  # Check if any of the expected warning patterns are present
  expect_true(
    any(grepl("Providing the sample size", warnings)) || 
    any(grepl("did not coverage", warnings))
  )
  
  # Still should get a valid result
  expect_s3_class(result, "colocboost")
})

# Test 8: Handling of NA variants
test_that("colocboost removes NA variants correctly", {
  sumstat_with_na <- test_sumstat_data$sumstat
  sumstat_with_na[[1]]$variant[1:2] <- NA
  
  warnings <- capture_warnings({
    validated_data <- colocboost_validate_input_data(
      sumstat = sumstat_with_na,
      LD = test_sumstat_data$LD
    )
  })
  
  expect_true(any(grepl("Removed variant with NA from sumstat 1", warnings)))
  # Should have 2 fewer variants
  expect_true(length(validated_data$Z[[1]]) == length(validated_data$Z[[2]])-2)
})

# Test 9: Handling of duplicate variants
test_that("colocboost removes duplicate variants correctly", {
  sumstat_with_dup <- test_sumstat_data$sumstat
  # Duplicate the first variant
  sumstat_with_dup[[1]] <- rbind(
    sumstat_with_dup[[1]][1, ],
    sumstat_with_dup[[1]]
  )
  
  warnings <- capture_warnings({
    result <- colocboost(
      sumstat = sumstat_with_dup,
      LD = test_sumstat_data$LD
    )
  })
  
  expect_true(any(grepl("Removed duplicate variants from sumstat 1", warnings)))
  expect_s3_class(result, "colocboost")
  # Should be back to original count
  expect_equal(length(result$data_info$variables), ncol(test_sumstat_data$X))
})

# Test 10: Handling of mismatched variants between sumstat and LD
test_that("colocboost handles variant mismatch between sumstat and LD", {
  # Create LD with different variant names
  LD_modified <- test_sumstat_data$LD
  rownames(LD_modified) <- colnames(LD_modified) <- paste0("SNP", 11:30)
  
  warnings <- capture_warnings({
    result <- colocboost_validate_input_data(
      sumstat = test_sumstat_data$sumstat,
      LD = LD_modified
    )
  })
  
  # Should warn about removing variants
  expect_true(any(grepl("removing.*variants since those variants are not in LD", warnings)))
})

# Test 11: Error when no common variants
test_that("colocboost errors with no common variants", {
  # Create LD with completely different variant names
  LD_no_overlap <- test_sumstat_data$LD
  rownames(LD_no_overlap) <- colnames(LD_no_overlap) <- paste0("DIFF", 1:20)
  
  warnings <- capture_warnings(
    colocboost_validate_input_data(
      sumstat = test_sumstat_data$sumstat,
      LD = LD_no_overlap
    )
  )
  expect_true(any(grepl("removing.*variants since those variants are not in LD", warnings)))
  expect_true(any(grepl("is empty after filtering", warnings)))
})

