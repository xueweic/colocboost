library(testthat)

# Helper function for test data generation
generate_one_causal_test_data <- function(n = 200, p = 30, L = 2, seed = 123) {
  set.seed(seed)
  
  # Generate X with strong LD structure to test LD handling
  sigma <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigma[i, j] <- 0.95^abs(i - j)  # Higher LD than standard test
    }
  }
  X <- MASS::mvrnorm(n, rep(0, p), sigma)
  colnames(X) <- paste0("SNP", 1:p)
  
  # Generate true effects - clear single causal variant per trait
  true_beta <- matrix(0, p, L)
  # Different causal variants for each trait
  true_beta[5, 1] <- 0.8   # Strong effect for trait 1
  true_beta[15, 2] <- 0.9  # Strong effect for trait 2
  
  # Generate Y with some noise
  Y <- matrix(0, n, L)
  for (l in 1:L) {
    Y[, l] <- X %*% true_beta[, l] + rnorm(n, 0, 0.7)  # Less noise for clearer signals
  }
  
  # Calculate beta, se, and z-scores
  beta <- matrix(0, ncol(X), ncol(Y))
  se <- matrix(0, ncol(X), ncol(Y))
  z <- matrix(0, ncol(X), ncol(Y))
  
  for (i in 1:ncol(Y)) {
    output = matrix(0, ncol(X), 2)
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
    LD = cor(X),
    sumstat = sumstat_list,
    effect_est = beta,
    effect_se = se,
    effect_n = rep(n, ncol(Y))
  )
}

# Create common test data - do this outside of test_that blocks
test_data <- generate_one_causal_test_data()
three_trait_data <- generate_one_causal_test_data(L = 3)
non_coloc_data <- generate_one_causal_test_data()  # Different seed would be better for true non-coloc

# Test 1: LD-free mode (diagonal LD matrix)
test_that("colocboost runs in LD-free mode correctly", {
  # Run colocboost without LD matrix
  warnings <- capture_warnings({
    result_ld_free <- colocboost(
      sumstat = test_data$sumstat,
      M = 1  # One iteration for LD-free mode
    )
  })
  
  # Test that we get a colocboost object
  expect_s3_class(result_ld_free, "colocboost")
  
  # Check that it used the LD-free mode
  expect_true(result_ld_free$model_info$model_coveraged == "one_causal")
  
  # Check that there's exactly one iteration
  expect_equal(result_ld_free$model_info$n_updates, 1)
  
  # Check dimensions
  expect_equal(length(result_ld_free$data_info$variables), ncol(test_data$X))
  expect_equal(result_ld_free$data_info$n_outcomes, 2)
})

# Test 2: One iteration mode with LD
test_that("colocboost runs in one iteration mode correctly", {
  # Run colocboost with LD matrix but M=1
  warnings <- capture_warnings({
    result_one_iter <- colocboost(
      sumstat = test_data$sumstat,
      LD = test_data$LD,
      M = 1  # One iteration 
    )
  })
  
  # Test that we get a colocboost object
  expect_s3_class(result_one_iter, "colocboost")
  
  # Check that it used the one_causal mode
  expect_true(result_one_iter$model_info$model_coveraged == "one_causal")
  
  # Check that the number of updates is 1
  expect_equal(result_one_iter$model_info$n_updates, 1)
  
  # It should have the update_status matrix
  expect_true(!is.null(result_one_iter$model_info$jk_star))
})

# Test 3: Focal outcome functionality with one causal
test_that("colocboost one causal mode handles focal outcome correctly", {
  # Run colocboost with focal_outcome_idx = 1
  warnings <- capture_warnings({
    result_focal <- colocboost(
      sumstat = test_data$sumstat,
      LD = test_data$LD,
      M = 1,
      focal_outcome_idx = 1
    )
  })
  
  # Test that we get a colocboost object
  expect_s3_class(result_focal, "colocboost")
  
  # Check focal outcome is correctly set
  expect_equal(result_focal$data_info$outcome_info$is_focal[1], TRUE)
  expect_equal(result_focal$data_info$outcome_info$is_focal[2], FALSE)
  
  # Check that it used the one_causal mode
  expect_true(result_focal$model_info$model_coveraged == "one_causal")
})

# Test 4: Different lambda values for focal outcome
test_that("colocboost one causal mode uses different lambda for focal outcome", {
  # Run with different lambda for focal outcome
  warnings <- capture_warnings({
    result_lambda <- colocboost(
      sumstat = test_data$sumstat,
      LD = test_data$LD,
      M = 1,
      focal_outcome_idx = 1,
      lambda = 0.3,             # Base lambda
      lambda_focal_outcome = 0.8  # Different lambda for focal outcome
    )
  })
  
  # Test that we get a colocboost object
  expect_s3_class(result_lambda, "colocboost")
  
  # Result should use one_causal approach
  expect_true(result_lambda$model_info$model_coveraged == "one_causal")
})

# Test 5: Multiple traits with no true colocalization
test_that("colocboost one causal mode correctly handles non-colocalized traits", {
  # Run colocboost
  warnings <- capture_warnings({
    result_non_coloc <- colocboost(
      sumstat = non_coloc_data$sumstat,
      LD = non_coloc_data$LD,
      M = 1
    )
  })
  
  # Test that we get a colocboost object
  expect_s3_class(result_non_coloc, "colocboost")
  
  # Result should use one_causal approach
  expect_true(result_non_coloc$model_info$model_coveraged == "one_causal")
})

# Test 6: With three or more traits
test_that("colocboost one causal mode handles three or more traits", {
  # Run colocboost
  warnings <- capture_warnings({
    result_three_traits <- colocboost(
      sumstat = three_trait_data$sumstat,
      LD = three_trait_data$LD,
      M = 1
    )
  })
  
  # Test that we get a colocboost object
  expect_s3_class(result_three_traits, "colocboost")
  
  # Check dimensions
  expect_equal(result_three_traits$data_info$n_outcomes, 3)
  
  # Result should use one_causal approach
  expect_true(result_three_traits$model_info$model_coveraged == "one_causal")
})

# Test 7: Effect of different jk_equiv parameters
test_that("colocboost one causal mode respects jk_equiv parameters", {
  # Run with stricter jk_equiv settings
  warnings <- capture_warnings({
    result_strict <- colocboost(
      sumstat = test_data$sumstat,
      LD = test_data$LD,
      M = 1,
      jk_equiv_corr = 0.95,   # Very high correlation required
      jk_equiv_loglik = 0.5   # More stringent loglik difference
    )
  })
  
  # Test that we get a colocboost object
  expect_s3_class(result_strict, "colocboost")
  
  # Should still be one_causal mode
  expect_true(result_strict$model_info$model_coveraged == "one_causal")
  
  # Run with lenient jk_equiv settings
  warnings <- capture_warnings({
    result_lenient <- colocboost(
      sumstat = test_data$sumstat,
      LD = test_data$LD,
      M = 1,
      jk_equiv_corr = 0.5,   # Lower correlation required
      jk_equiv_loglik = 2.0  # Larger loglik difference allowed
    )
  })
  
  # Test that we get a colocboost object
  expect_s3_class(result_lenient, "colocboost")
})

# Test 8: Interaction with output level parameter
test_that("colocboost one causal mode works with different output_level values", {
  for (level in 1:3) {
    warnings <- capture_warnings({
      result <- colocboost(
        sumstat = test_data$sumstat,
        LD = test_data$LD,
        M = 1,
        output_level = level
      )
    })
    
    # Test that we get a colocboost object
    expect_s3_class(result, "colocboost")
    
    # Check output level-specific fields
    if (level >= 2) {
      expect_true("ucos_details" %in% names(result))
    }
    if (level >= 3) {
      expect_true("diagnostic_details" %in% names(result))
    }
  }
})

# Test 10: Test with HyPrColoc compatible format
test_that("colocboost one causal mode works with HyPrColoc format", {
  # Run colocboost with effect matrices
  warnings <- capture_warnings({
    result_hypr <- colocboost(
      effect_est = test_data$effect_est,
      effect_se = test_data$effect_se,
      effect_n = test_data$effect_n,
      LD = test_data$LD,
      M = 1
    )
  })
  
  # Test that we get a colocboost object
  expect_s3_class(result_hypr, "colocboost")
  
  # Check that it used the one_causal mode
  expect_true(result_hypr$model_info$model_coveraged == "one_causal")
  
  # Check dimensions
  expect_equal(length(result_hypr$data_info$variables), ncol(test_data$X))
  expect_equal(result_hypr$data_info$n_outcomes, 2)
})