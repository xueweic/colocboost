library(testthat)

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

# Test for check_jk_jkeach function
test_that("check_jk_jkeach identifies equivalent variants based on LD", {
  skip("Internal function being tested through main function")
  
  # This function is harder to test directly as it relies on internal colocboost objects
  # We'll test it indirectly through colocboost main function
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

# Test for get_max_profile function
test_that("get_max_profile updates check_null_max correctly", {
  skip("Internal function being tested through main function")
  
  # This function modifies internal CB object - hard to test directly
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
