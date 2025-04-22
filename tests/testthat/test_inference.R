library(testthat)

# Utility function to generate a simple colocboost results 
generate_test_result <- function(n = 100, p = 50, L = 2, seed = 42) {
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
    true_beta[40, 1] <- 0.6 # SNP10 also affects the trait
  } else {
    # Multi-trait case
    true_beta[5, 1] <- 0.7  # SNP5 affects trait 1
    true_beta[5, 2] <- 0.6  # SNP5 also affects trait 2 (colocalized)
    true_beta[40, 2] <- 0.8 # SNP10 only affects trait 2
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

# Test for get_strong_colocalization
test_that("get_robust_colocalization filters results correctly", {

  # Generate a test colocboost results
  cb_res <- generate_test_result()
  
  # Basic call
  expect_error(get_robust_colocalization(cb_res), NA)
  
  # With stricter thresholds
  expect_error(get_robust_colocalization(cb_res, cos_npc_cutoff = 0.8), NA)
  
  # With p-value threshold
  expect_error(get_robust_colocalization(cb_res, pvalue_cutoff = 0.05), NA)
})

# Test for get_hierarchical_clusters
test_that("get_hierarchical_clusters functions correctly", {
  # Test case 1: Simple 2x2 correlation matrix with high correlation
  set.seed(123)
  cormat_high <- matrix(c(1, 0.95, 0.95, 1), nrow = 2)
  result_high <- get_hierarchical_clusters(cormat_high)
  
  # Should return a single cluster since correlation > min_cluster_corr
  expect_equal(ncol(result_high$cluster), 1)
  expect_equal(nrow(result_high$cluster), 2)
  expect_equal(sum(result_high$cluster), 2) # All items in one cluster
  
  # Test case 2: Matrix with two distinct clusters
  set.seed(456)
  N <- 100
  P <- 4
  # Generate correlation structure with two distinct clusters
  sigma <- matrix(0.2, nrow = P, ncol = P)
  diag(sigma) <- 1
  sigma[1:2, 1:2] <- 0.9
  sigma[3:4, 3:4] <- 0.9
  X <- MASS::mvrnorm(N, rep(0, P), sigma)
  cormat_two <- cor(X)
  
  result_two <- get_hierarchical_clusters(cormat_two)
  
  # Should identify two clusters
  expect_equal(ncol(result_two$cluster), 2)
  expect_equal(nrow(result_two$cluster), P)
  
  # Test case 3: Single variable edge case
  cormat_single <- matrix(1, nrow = 1)
  result_single <- get_hierarchical_clusters(cormat_single)
  
  expect_equal(ncol(result_single$cluster), 1)
  expect_equal(nrow(result_single$cluster), 1)
  
  # Test case 4: Larger matrix with multiple clusters
  set.seed(789)
  P <- 10
  sigma_large <- matrix(0.1, nrow = P, ncol = P)
  diag(sigma_large) <- 1
  # Create 3 distinct clusters
  sigma_large[1:3, 1:3] <- 0.8
  sigma_large[4:6, 4:6] <- 0.8
  sigma_large[7:10, 7:10] <- 0.8
  X_large <- MASS::mvrnorm(N, rep(0, P), sigma_large)
  cormat_large <- cor(X_large)
  
  result_large <- get_hierarchical_clusters(cormat_large, min_cluster_corr = 0.7)
  
  # Should identify 3 clusters
  expect_equal(ncol(result_large$cluster), 3)
  expect_equal(nrow(result_large$cluster), P)
  
  # Test case 5: Test with different min_cluster_corr
  result_diff_threshold <- get_hierarchical_clusters(cormat_two, min_cluster_corr = 0.5)
  
  # With lower threshold, might merge clusters
  expect_true(ncol(result_diff_threshold$cluster) <= 2)
  
  # Test get_modularity directly
  # Test case 6: Test modularity calculation
  B <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1), ncol = 2)
  W <- matrix(c(
    1.0, 0.9, 0.2, 0.1,
    0.9, 1.0, 0.3, 0.2,
    0.2, 0.3, 1.0, 0.8,
    0.1, 0.2, 0.8, 1.0
  ), nrow = 4)
  
  modularity <- get_modularity(W, B)
  expect_type(modularity, "double")
  expect_true(!is.na(modularity))
  
  # Test case 7: Edge cases for get_modularity
  # Empty matrix
  W_empty <- matrix(0, nrow = 2, ncol = 2)
  B_empty <- matrix(c(1, 1), ncol = 1)
  mod_empty <- get_modularity(W_empty, B_empty)
  expect_equal(mod_empty, 0)
  
  # Matrix with only positive weights
  W_pos <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  B_pos <- matrix(c(1, 1), ncol = 1)
  mod_pos <- get_modularity(W_pos, B_pos)
  expect_true(!is.na(mod_pos))
  
  # Matrix with only negative weights
  W_neg <- matrix(c(1, -0.5, -0.5, 1), nrow = 2)
  B_neg <- matrix(c(1, 1), ncol = 1)
  mod_neg <- get_modularity(W_neg, B_neg)
  expect_true(!is.na(mod_neg))
  
  # Test get_n_cluster
  # Test case 8: Test cluster number determination
  hc <- hclust(as.dist(1 - cormat_two))
  cluster_result <- get_n_cluster(hc, cormat_two)
  
  expect_type(cluster_result, "list")
  expect_true("n_cluster" %in% names(cluster_result))
  expect_true("Qmodularity" %in% names(cluster_result))
  
  # Test case 9: All correlations above threshold
  cormat_all_high <- matrix(0.9, nrow = 3, ncol = 3)
  diag(cormat_all_high) <- 1
  hc_high <- hclust(as.dist(1 - cormat_all_high))
  result_all_high <- get_n_cluster(hc_high, cormat_all_high, min_cluster_corr = 0.85)
  
  expect_equal(result_all_high$n_cluster, 1)
})


# Test get_robust_colocalization function
test_that("get_ambiguous_colocalization identifies ambiguous colocalizations correctly", {
  # The function expects a specialized test dataset that has ambiguous colocalizations
  # There's a reference in the example to a dataset named "Ambiguous_Colocalization"
  data(Ambiguous_Colocalization)
  test_colocboost_results <- Ambiguous_Colocalization$ColocBoost_Results
  
  # Basic call with default parameters
  result <- get_ambiguous_colocalization(test_colocboost_results)
  
  # Check that the returned object is of class "colocboost"
  expect_s3_class(result, "colocboost")
  
  # Check that the ambiguous_ucos field exists in the result
  expect_true("ambigous_ucos" %in% names(result))
  
  # If ambiguous colocalizations were found, test their structure
  if (length(result$ambigous_ucos) > 0) {
    # There should be fields for the ambiguous UCOs details
    expect_true("ambigouse_ucos" %in% names(result$ambigous_ucos[[1]]))
    expect_true("ambigouse_ucos_outcomes" %in% names(result$ambigous_ucos[[1]]))
    expect_true("ambigous_ucos_weight" %in% names(result$ambigous_ucos[[1]]))
    expect_true("ambigous_ucos_puriry" %in% names(result$ambigous_ucos[[1]]))
    expect_true("ambigouse_ucos_unoin" %in% names(result$ambigous_ucos[[1]]))
    expect_true("ambigouse_ucos_overlap" %in% names(result$ambigous_ucos[[1]]))
  }
  
  # Test with custom correlation thresholds
  result_custom <- get_ambiguous_colocalization(
    test_colocboost_results,
    min_abs_corr_between_ucos = 0.7,
    median_abs_corr_between_ucos = 0.9
  )
  
  # The result should still be a colocboost object
  expect_s3_class(result_custom, "colocboost")
  
  # Test with input that has no ucos_details
  # Create a modified object without ucos_details
  test_no_ucos <- test_colocboost_results
  test_no_ucos$ucos_details <- NULL
  
  # Should show a warning but not error
  expect_warning(get_ambiguous_colocalization(test_no_ucos))
  
  # Test with input that has only one trait-specific effect
  cb_res <- generate_test_result(n=500)
  
  expect_message(
    result <- get_ambiguous_colocalization(cb_res),
    "Only one trait-specific \\(uncolocalized\\) effect in this region!"
  )
  
  # The result should be unchanged from the input
  expect_equal(result, cb_res)
  
  # There should be no ambiguous_ucos field added
  expect_false("ambigous_ucos" %in% names(result))

})

