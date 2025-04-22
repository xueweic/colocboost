library(testthat)

# Utility function to generate a simple colocboost results 
generate_test_result <- function(n = 100, p = 20, L = 2, seed = 42) {
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
  
  # Generate true effects based on the number of traits
  true_beta <- matrix(0, p, L)
  
  if (L == 1) {
    # Single trait case
    true_beta[5, 1] <- 0.7  # SNP5 affects the trait
    true_beta[10, 1] <- 0.6 # SNP10 also affects the trait
  } else {
    # Multi-trait case
    true_beta[5, 1] <- 0.7  # SNP5 affects trait 1
    true_beta[5, 2] <- 0.6  # SNP5 also affects trait 2 (colocalized)
    true_beta[10, 2] <- 0.5 # SNP10 only affects trait 2
  }
  
  # Generate Y with some noise
  Y <- matrix(0, n, L)
  for (l in 1:L) {
    Y[, l] <- X %*% true_beta[, l] + rnorm(n, 0, 1)
  }
  
  # Prepare input for colocboost based on L
  if (L == 1) {
    # For single trait, Y should be a vector
    Y_input <- Y[,1]
    X_input <- X
  } else {
    # For multiple traits, convert to list format
    Y_input <- lapply(1:L, function(l) Y[,l])
    X_input <- replicate(L, X, simplify = FALSE)
  }
  
  # Run colocboost with minimal parameters to get a model object
  suppressWarnings({
    result <- colocboost(
      X = X_input, 
      Y = Y_input,
      M = 5,  # Small number of iterations for faster testing
      output_level = 3  # Include full model details
    )
  })
  
  return(result)
}

# Test for get_integrated_weight function
test_that("get_integrated_weight correctly integrates weights", {
  # Create test weights matrix
  weights <- matrix(c(0.3, 0.2, 0.1, 0.4, 0.4, 0.3, 0.2, 0.1), nrow = 4)
  colnames(weights) <- c("outcome1", "outcome2")
  
  # Call the function
  result <- get_integrated_weight(weights, weight_fudge_factor = 1.5)
  
  # Test expected properties
  expect_equal(length(result), nrow(weights))
  expect_equal(sum(result), 1, tolerance = 1e-10)  # Should sum to 1
  
  # Higher values in both columns should get higher integrated weights
  expect_true(result[1] > result[3])
  expect_true(result[4] < result[2])
})

# Test for w_cs function
test_that("w_cs correctly identifies confidence set for weight vector", {
  # Create test weight vector
  w <- c(0.5, 0.3, 0.1, 0.1)
  
  # Call function
  result <- w_cs(w, coverage = 0.8)
  
  # Expected result for 80% coverage: w[1] + w[2] = 0.8, so first 2 elements should be 1
  expected <- c(1, 1, 0, 0)
  expect_equal(result, expected)
  
  # Test with different coverage
  result2 <- w_cs(w, coverage = 0.9)
  expected2 <- c(1, 1, 1, 1)  # First 4 elements cover 90% since 3 and 4 with the same 0.1 weight 
  expect_equal(result2, expected2)
})

# Test for get_in_cos function
test_that("get_in_cos correctly identifies indices in confidence set", {
  # Create test weight vector
  w <- c(0.1, 0.5, 0.2, 0.1, 0.1)
  
  # Call function
  result <- get_in_cos(w, coverage = 0.7)
  
  # Expected result: top elements that sum to 0.7 coverage
  # Ordering: 2, 3, 1, 4, 5
  # So elements 2 and 3 (indexes 2 and 3) should be in the CS
  expect_equal(result[[1]], c(2, 3))
  
  # Test with higher coverage
  result2 <- get_in_cos(w, coverage = 0.9)
  # Should include elements 2, 3, 1, 4, 5
  expect_equal(sort(result2[[1]]), c(1, 2, 3, 4, 5))
})

# Test for merge_sets function
test_that("merge_sets correctly merges overlapping sets", {
  # Create test vector of semicolon-delimited sets
  vec <- c("1;2;3", "3;4;5", "6;7", "8;9", "9;10", "11")
  
  # Call function
  result <- merge_sets(vec)
  
  # Expected result: Sets {1,2,3,4,5}, {6,7}, {8,9,10}, {11}
  expect_equal(length(result), 4)
  expect_true("1;2;3;4;5" %in% result)
  expect_true("6;7" %in% result)
  expect_true("8;9;10" %in% result)
  expect_true("11" %in% result)
})

# Test for get_purity function
test_that("get_purity calculates correlation statistics correctly", {
  # Create test data with known correlation structure
  set.seed(123)
  X <- matrix(rnorm(100*3), 100, 3)
  X[,2] <- X[,1] + rnorm(100, 0, 0.2)  # Column 2 highly correlated with column 1
  X[,3] <- rnorm(100)  # Column 3 independent
  
  Xcorr <- cor(X)
  
  # Test with a single position
  result_single <- get_purity(1, X = X)
  expect_equal(result_single, c(1, 1, 1))  # For single element, all correlations are 1
  
  # Test with highly correlated positions
  result_corr <- get_purity(c(1, 2), X = X)
  expect_gt(result_corr[1], 0.8)  # Minimum correlation should be high
  
  # Test with uncorrelated positions
  result_uncorr <- get_purity(c(1, 3), X = X)
  expect_lt(result_uncorr[1], 0.5)  # Minimum correlation should be low
  
  # Test with XCorr input
  result_Xcorr <- get_purity(c(1, 2), Xcorr = Xcorr)
  expect_gt(result_Xcorr[1], 0.8)  # Should give similar results to X input
})

# Test for get_cormat function
test_that("get_cormat calculates correlation matrix correctly", {
  # Create test data
  set.seed(123)
  X <- matrix(rnorm(100*5), 100, 5)
  
  # Make some columns correlated
  X[,2] <- X[,1] + rnorm(100, 0, 0.2)
  
  # Calculate correlation matrix with base R
  expected <- cor(X)
  
  # Calculate with get_cormat
  result <- get_cormat(X)
  
  # Should match the R correlation matrix
  expect_equal(result, expected, tolerance = 1e-6)
})

# Test for get_between_purity function
test_that("get_between_purity calculates correlation between sets", {
  # Create test data with known correlation structure
  set.seed(123)
  X <- matrix(rnorm(100*5), 100, 5)
  X[,2] <- X[,1] + rnorm(100, 0, 0.2)  # Column 2 highly correlated with column 1
  X[,4] <- X[,3] + rnorm(100, 0, 0.2)  # Column 4 highly correlated with column 3
  
  Xcorr <- cor(X)
  
  # Test between correlated sets
  set1 <- c(1, 2)
  set2 <- c(1, 3)
  result <- get_between_purity(set1, set2, X = X)
  
  # Min, max, median correlations
  expect_length(result, 3)
  
  # Check that correlations between sets with shared elements are high
  expect_gt(result[2], 0.8)  # Max correlation should be high (same element)
  
  # Test between uncorrelated sets
  set1 <- c(1, 2)
  set2 <- c(3, 5)
  result2 <- get_between_purity(set1, set2, X = X)
  
  # Min correlation should be low
  expect_lt(result2[1], 0.5)
})


# Test for get_merge_ordered_with_indices function
test_that("get_merge_ordered_with_indices merges vectors", {
  # Create test vectors
  vec1 <- c("a", "b", "c")
  vec2 <- c("b", "c", "d")
  vec3 <- c("c", "d", "e")
  
  vector_list <- list(vec1, vec2, vec3)
  
  # Call function
  result <- get_merge_ordered_with_indices(vector_list)
  
  # Expected result: a, b, c, d, e (unique elements in order of appearance)
  expect_equal(result, c("a", "b", "c", "d", "e"))
  
  # Test with duplicate vectors
  vector_list2 <- list(vec1, vec1, vec2, vec3)
  result2 <- get_merge_ordered_with_indices(vector_list2)
  expect_equal(result2, c("a", "b", "c", "d", "e"))
})

# Test for get_cos_purity function
test_that("get_cos_purity calculates correct purity statistics", {
  # Set up test data
  set.seed(123)
  N <- 100
  P <- 20
  
  # Generate X with correlation structure
  sigma <- 0.8^abs(outer(1:P, 1:P, "-"))
  X <- MASS::mvrnorm(N, rep(0, P), sigma)
  colnames(X) <- paste0("SNP", 1:P)
  
  # Create test CoS lists
  cos1 <- c(1, 2, 3)
  cos2 <- c(10, 11, 12)
  cos3 <- c(5, 6, 7, 8)
  cos_list <- list(trait1 = cos1, trait2 = cos2, trait3 = cos3)
  
  # Calculate X correlation matrix for testing both methods
  Xcorr <- cor(X)
  
  # Test with X matrix input
  result_X <- get_cos_purity(cos_list, X = X)
  
  # Test with correlation matrix input
  result_Xcorr <- get_cos_purity(cos_list, Xcorr = Xcorr)
  
  # Basic structure tests
  expect_type(result_X, "list")
  expect_equal(length(result_X), 3)
  expect_named(result_X, c("min_abs_cor", "max_abs_cor", "median_abs_cor"))
  
  # Check dimensions
  expect_equal(dim(result_X$min_abs_cor), c(3, 3))
  expect_equal(rownames(result_X$min_abs_cor), c("trait1", "trait2", "trait3"))
  
  # Test symmetry of results
  expect_equal(result_X$min_abs_cor[1, 2], result_X$min_abs_cor[2, 1])
  expect_equal(result_X$max_abs_cor[1, 3], result_X$max_abs_cor[3, 1])
  
  # Test that X and Xcorr methods give same results (within numerical precision)
  expect_equal(result_X$min_abs_cor, result_Xcorr$min_abs_cor, tolerance = 1e-10)
  expect_equal(result_X$max_abs_cor, result_Xcorr$max_abs_cor, tolerance = 1e-10)
  expect_equal(result_X$median_abs_cor, result_Xcorr$median_abs_cor, tolerance = 1e-10)
  
  # Test single CoS case
  single_cos <- list(single = cos1)
  result_single <- get_cos_purity(single_cos, X = X)
  expect_equal(dim(result_single$min_abs_cor), c(1, 1))
  expect_named(result_single, c("min_abs_cor", "max_abs_cor", "median_abs_cor"))
  
  # Test empty CoS case
  expect_null(get_cos_purity(list(), X = X))
  
  # Test n_purity parameter
  result_limited <- get_cos_purity(cos_list, X = X, n_purity = 2)
  # The result still has the same structure
  expect_type(result_limited, "list")
  expect_equal(length(result_limited), 3)
  
  # Test error cases
  expect_error(get_cos_purity(cos_list), "Either X or Xcorr must be provided")
  
  # Test converting numeric input to list
  num_cos <- 1:5
  result_num <- get_cos_purity(num_cos, X = X)
  expect_equal(dim(result_num$min_abs_cor), c(1, 1))
})

# Test for get_cos function
test_that("get_cos extracts CoS correctly with generated test results", {
  # Generate test data
  cb_output <- generate_test_result()
  
  # Test basic functionality with default coverage
  result_default <- get_cos(cb_output, coverage = 0.95)
  
  # Check structure
  expect_type(result_default, "list")
  expect_named(result_default, c("cos", "cos_purity"))
  expect_type(result_default$cos, "list")
  expect_null(result_default$cos_purity)  # Should be NULL since X and Xcorr not provided
  
  # Get names of CoS in the result
  cos_names <- names(result_default$cos)
  expect_true(all(grepl("_coverage_0.95$", cos_names)))  # Names should include coverage value
  
  # Test with different coverage values
  result_high <- get_cos(cb_output, coverage = 0.99)
  result_low <- get_cos(cb_output, coverage = 0.75)
  
  # Higher coverage should include more SNPs
  for (cos_name_base in gsub("_coverage_0.95$", "", cos_names)) {
    high_cos_name <- paste0(cos_name_base, "_coverage_0.99")
    low_cos_name <- paste0(cos_name_base, "_coverage_0.75")
    default_cos_name <- paste0(cos_name_base, "_coverage_0.95")
    
    # Check if names exist in all results
    if (high_cos_name %in% names(result_high$cos) && 
        low_cos_name %in% names(result_low$cos)) {
      # Higher coverage should have equal or more SNPs
      expect_true(
        length(result_high$cos[[high_cos_name]]) >= 
        length(result_default$cos[[default_cos_name]])
      )
      # Lower coverage should have equal or fewer SNPs
      expect_true(
        length(result_low$cos[[low_cos_name]]) <= 
        length(result_default$cos[[default_cos_name]])
      )
    }
  }
  
  # Create X and Xcorr for further testing
  # Extract SNP names from the first CoS
  first_cos <- result_default$cos[[1]]
  all_snps <- unique(unlist(result_default$cos))
  max_snp_index <- max(as.integer(gsub("SNP", "", all_snps)))
  
  # Create a small genotype matrix with the necessary SNPs
  set.seed(123)
  N <- 100
  P <- max(max_snp_index + 10, 50)  # Ensure P is large enough
  
  # Generate X with LD structure
  sigma <- 0.8^abs(outer(1:P, 1:P, "-"))
  X <- MASS::mvrnorm(N, rep(0, P), sigma)
  colnames(X) <- paste0("SNP", 1:P)
  
  # Calculate correlation matrix
  Xcorr <- cor(X)
  
  # Test with X provided
  result_with_X <- get_cos(cb_output, coverage = 0.95, X = X)
  
  # Check structure with X provided
  expect_type(result_with_X, "list")
  expect_named(result_with_X, c("cos", "cos_purity"))
  expect_type(result_with_X$cos, "list")
  expect_type(result_with_X$cos_purity, "list")
  
  # Check purity matrices
  expect_named(result_with_X$cos_purity, c("min_abs_cor", "max_abs_cor", "median_abs_cor"))
  expect_true(all(sapply(result_with_X$cos_purity, is.matrix)))
  
  # Test with Xcorr provided
  result_with_Xcorr <- get_cos(cb_output, coverage = 0.95, Xcorr = Xcorr)
  
  # Check structure with Xcorr provided
  expect_type(result_with_Xcorr, "list")
  expect_named(result_with_Xcorr, c("cos", "cos_purity"))
  expect_type(result_with_Xcorr$cos, "list")
  expect_type(result_with_Xcorr$cos_purity, "list")
  
  # Test with min_abs_corr filtering
  result_high_purity <- get_cos(cb_output, coverage = 0.95, X = X, min_abs_corr = 0.8)
  
  # High purity threshold might remove some CoS
  expect_type(result_high_purity, "list")
  expect_named(result_high_purity, c("cos", "cos_purity"))
  
  # If result_high_purity$cos is not NULL, it should have fewer or equal CoS compared to result_with_X
  if (!is.null(result_high_purity$cos)) {
    expect_true(length(result_high_purity$cos) <= length(result_with_X$cos))
  }
  
  # Test with median_abs_corr filtering
  result_median_purity <- get_cos(cb_output, coverage = 0.95, X = X, min_abs_corr = 0.4, median_abs_corr = 0.7)
  
  # Should have structure similar to other results
  expect_type(result_median_purity, "list")
  expect_named(result_median_purity, c("cos", "cos_purity"))
  
  # Test empty colocalization results
  empty_cb_output <- cb_output
  empty_cb_output$cos_details$cos <- NULL
  
  expect_warning(
    result_empty <- get_cos(empty_cb_output, coverage = 0.95),
    "No colocalization results in this region!"
  )
  expect_equal(result_empty, list("cos" = NULL, "cos_purity" = NULL))
})