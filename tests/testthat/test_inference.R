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


# Utility function to generate test data with uncolocalized effects
generate_ucos_test_data <- function(n = 500, p = 60, L = 3, seed = 42, output_level = 3) {
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
  
  # Generate true effects with both colocalized and trait-specific effects
  true_beta <- matrix(0, p, L)
  
  # Colocalized effect: SNP10 affects traits 1 and 2
  true_beta[10, 1] <- 0.7
  true_beta[10, 2] <- 0.6
  
  # Trait-specific (uncolocalized) effects
  true_beta[30, 1] <- 0.5  # SNP30 only affects trait 1
  true_beta[45, 2] <- 0.6  # SNP45 only affects trait 2
  true_beta[50, 3] <- 0.7  # SNP50 only affects trait 3
  
  # Generate Y with some noise
  Y <- matrix(0, n, L)
  for (l in 1:L) {
    Y[, l] <- X %*% true_beta[, l] + rnorm(n, 0, 1)
  }
  
  # Prepare input for colocboost
  Y_input <- lapply(1:L, function(l) Y[,l])
  X_input <- replicate(L, X, simplify = FALSE)
  
  # Run colocboost with output_level to get ucos_details
  suppressWarnings({
    result <- colocboost(
      X = X_input, 
      Y = Y_input,
      output_level = output_level
    )
  })
  
  return(result)
}



# Test for get_strong_colocalization
test_that("get_robust_colocalization filters results correctly", {

  # Generate a test colocboost results
  cb_res <- generate_test_result()
  
  # Basic call
  expect_error(suppressWarnings(get_robust_colocalization(cb_res)), NA)
  
  # With stricter thresholds
  expect_error(suppressWarnings(get_robust_colocalization(cb_res, cos_npc_cutoff = 0.8)), NA)
  
  # With p-value threshold
  expect_error(suppressWarnings(get_robust_colocalization(cb_res, pvalue_cutoff = 0.05)), NA)
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


# Test get_ambiguous_colocalization function
test_that("get_ambiguous_colocalization identifies ambiguous colocalizations correctly", {
  # The function expects a specialized test dataset that has ambiguous colocalizations
  # There's a reference in the example to a dataset named "Ambiguous_Colocalization"
  data("Ambiguous_Colocalization", package = "colocboost", envir = environment())
  test_colocboost_results <- Ambiguous_Colocalization$ColocBoost_Results
  
  # Basic call with default parameters
  result <- get_ambiguous_colocalization(test_colocboost_results)
  
  # Check that the returned object is of class "colocboost"
  expect_s3_class(result, "colocboost")
  
  # Check that the ambiguous_cos field exists in the result
  expect_true("ambiguous_cos" %in% names(result))
  
  # If ambiguous colocalizations were found, test their structure
  if (length(result$ambigous_cos) > 0) {
    # There should be fields for the ambiguous UCOs details
    expect_true("ambiguous_cos" %in% names(result$ambigous_ucos[[1]]))
    expect_true("ambiguous_cos_overlap" %in% names(result$ambigous_ucos[[1]]))
    expect_true("ambiguous_cos_union" %in% names(result$ambigous_ucos[[1]]))
    expect_true("ambiguous_cos_outcomes" %in% names(result$ambigous_ucos[[1]]))
    expect_true("ambigous_cos_weight" %in% names(result$ambigous_ucos[[1]]))
    expect_true("ambigous_cos_purity" %in% names(result$ambigous_ucos[[1]]))
    expect_true("recalibrated_cos_vcp" %in% names(result$ambigous_ucos[[1]]))
    expect_true("recalibrated_cos" %in% names(result$ambigous_ucos[[1]]))
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
  
  # There should be no ambiguous_cos field added
  expect_false("ambigous_cos" %in% names(result))

})



# Test get_ucos_summary function
test_that("get_ucos_summary funtionality", {
  # The function expects a specialized test dataset that has ambiguous colocalizations
  # There's a reference in the example to a dataset named "Ambiguous_Colocalization"
  data("Ambiguous_Colocalization", package = "colocboost", envir = environment())
  test_colocboost_results <- Ambiguous_Colocalization$ColocBoost_Results
  
  # Basic call with default parameters
  summary <- get_ucos_summary(test_colocboost_results)
  
  # Check structure of summary table
  expect_true(is.data.frame(summary))
    
  # Check expected columns exist
  expected_cols <- c(
    "outcomes", "ucos_id", "purity",
    "top_variable", "top_variable_vpa", "n_variables", "ucos_index",
    "ucos_variables", "ucos_variables_vpa"
  )
  for (col in expected_cols) {
    expect_true(col %in% colnames(summary))
  }

  # Basic call with default parameters
  summary_ambiguous <- get_ucos_summary(test_colocboost_results, ambiguous_cos = TRUE)
  expect_true(all.equal(names(summary_ambiguous), c("ucos_summary", "ambiguous_cos_summary")))
    
  # Check expected columns exist
  expected_cols <- c(
    "outcomes", "ucos_id", "min_between_purity", "median_between_purity",
    "overlap_idx", "overlap_variables", "n_recalibrated_variables",
    "recalibrated_index", "recalibrated_variables", "recalibrated_variables_vcp"
  )
  for (col in expected_cols) {
    expect_true(col %in% colnames(summary_ambiguous$ambiguous_cos_summary))
  }

})


test_that("get_colocboost_summary works correctly", {
  # Setup mock data
  set.seed(1)
  N <- 1000
  P <- 100
  # Generate X with LD structure
  sigma <- 0.9^abs(outer(1:P, 1:P, "-"))
  X <- MASS::mvrnorm(N, rep(0, P), sigma)
  colnames(X) <- paste0("SNP", 1:P)
  L <- 3
  true_beta <- matrix(0, P, L)
  true_beta[10, 1] <- 0.5 # SNP10 affects trait 1
  true_beta[10, 2] <- 0.4 # SNP10 also affects trait 2 (colocalized)
  true_beta[50, 2] <- 0.3 # SNP50 only affects trait 2
  true_beta[80, 3] <- 0.6 # SNP80 only affects trait 3
  Y <- matrix(0, N, L)
  for (l in 1:L) {
    Y[, l] <- X %*% true_beta[, l] + rnorm(N, 0, 1)
  }
  
  # Run colocboost
  cb_output <- colocboost(X = X, Y = Y)
  
  # Test summary_level = 1 (default)
  summary1 <- get_colocboost_summary(cb_output)
  
  # Check structure
  expect_type(summary1, "list")
  expect_named(summary1, "cos_summary")
  expect_s3_class(summary1$cos_summary, "data.frame")
  
  # Check required columns in cos_summary
  expected_cols <- c("focal_outcome", "colocalized_outcomes", "cos_id", 
                    "purity", "top_variable", "top_variable_vcp", 
                    "cos_npc", "min_npc_outcome", "n_variables")
  expect_true(all(expected_cols %in% colnames(summary1$cos_summary)))
  
  # Test with outcome_names
  outcome_names <- c("Trait1", "Trait2", "Trait3")
  summary_with_names <- get_colocboost_summary(cb_output, outcome_names = outcome_names)
  coloc_outcome <- strsplit(summary_with_names$cos_summary$colocalized_outcomes, "; ")[[1]]
  expect_true(all( coloc_outcome %in% 
                 c("Trait1", "Trait2", "Trait3", paste(outcome_names, collapse = ", "))))
  
  # Test with region_name
  region_summary <- get_colocboost_summary(cb_output, region_name = "TestRegion")
  expect_true("region_name" %in% colnames(region_summary$cos_summary))
  expect_equal(unique(region_summary$cos_summary$region), "TestRegion")
  
  # Test summary_level = 2
  expect_warning(summary2 <- get_colocboost_summary(cb_output, summary_level = 2),
  "Please run colocboost model with output_level=2")

  cb_output <- colocboost(X = X, Y = Y, output_level = 2)
  summary2 <- get_colocboost_summary(cb_output, summary_level = 2)
  expect_named(summary2, c("cos_summary", "ucos_summary"))
  expect_s3_class(summary2$ucos_summary, "data.frame")
  
  # Test summary_level = 3
  summary3 <- get_colocboost_summary(cb_output, 
                                    summary_level = 3, 
                                    min_abs_corr_between_ucos = 0.4,
                                    median_abs_corr_between_ucos = 0.7)
  expect_named(summary3, c("cos_summary", "ucos_summary", "ambiguous_cos_summary"))
  expect_s3_class(summary3$ucos_summary, "data.frame")
  
  # Test with interest_outcome
  interest_summary <- get_colocboost_summary(cb_output, 
                                           outcome_names = outcome_names,
                                           interest_outcome = c("Trait1"))
  # Should only contain colocalization events involving Trait1
  if(nrow(interest_summary$cos_summary) > 0) {
    expect_true(all(sapply(interest_summary$cos_summary$colocalized_outcomes, 
                         function(x) grepl("Trait1", x))))
  }
  
  # Test error handling
  expect_error(get_colocboost_summary("not_a_colocboost_object"), 
              "Input must from colocboost output!")
})



# ============================================================================
# Tests for get_robust_ucos
# ============================================================================

test_that("get_robust_ucos basic functionality works", {
  
  # Generate test data with ucos
  cb_res <- generate_ucos_test_data(output_level = 2)
  
  # Skip if no ucos were detected
  skip_if(is.null(cb_res$ucos_details), "No ucos detected in test data")
  
  # Basic call with default parameters
  result <- get_robust_ucos(cb_res)
  
  # Should return a colocboost object
  expect_s3_class(result, "colocboost")
  
  # Should have ucos_details (unless all were filtered out)
  expect_true("ucos_details" %in% names(result))
})


test_that("get_robust_ucos filters by npc_outcome_cutoff", {
  
  # Generate test data
  cb_res <- generate_ucos_test_data(output_level = 2)
  
  # Skip if no ucos were detected
  skip_if(is.null(cb_res$ucos_details), "No ucos detected in test data")
  
  # Get original number of ucos
  n_ucos_original <- length(cb_res$ucos_details$ucos$ucos_index)
  
  # Apply lenient filtering (should keep most/all)
  result_lenient <- get_robust_ucos(cb_res, npc_outcome_cutoff = 0.1)
  
  # Apply strict filtering (should remove some/all)
  result_strict <- get_robust_ucos(cb_res, npc_outcome_cutoff = 0.8)
  
  # Check that strict filtering removes at least as many as lenient
  n_ucos_lenient <- if(is.null(result_lenient$ucos_details)) 0 else length(result_lenient$ucos_details$ucos$ucos_index)
  n_ucos_strict <- if(is.null(result_strict$ucos_details)) 0 else length(result_strict$ucos_details$ucos$ucos_index)
  
  expect_true(n_ucos_strict <= n_ucos_lenient)
  expect_true(n_ucos_lenient <= n_ucos_original)
})

test_that("get_robust_ucos filters by pvalue_cutoff", {
  
  # Generate test data
  cb_res <- generate_ucos_test_data(output_level = 2)
  
  # Skip if no ucos were detected
  skip_if(is.null(cb_res$ucos_details), "No ucos detected in test data")
  
  # Apply lenient p-value filtering
  expect_message(
    result_lenient <- get_robust_ucos(cb_res, npc_outcome_cutoff = 0, pvalue_cutoff = 0.1),
    "Extracting outcome-specific.*pvalue_cutoff = 0.1"
  )
  
  # Apply strict p-value filtering
  expect_message(
    result_strict <- get_robust_ucos(cb_res, npc_outcome_cutoff = 0, pvalue_cutoff = 1e-6),
    "Extracting outcome-specific.*pvalue_cutoff = 1e-06"
  )
  
  # Both should return colocboost objects
  expect_s3_class(result_lenient, "colocboost")
  expect_s3_class(result_strict, "colocboost")
})

test_that("get_robust_ucos filters by both npc_outcome_cutoff and pvalue_cutoff", {
  
  # Generate test data
  cb_res <- generate_ucos_test_data(output_level = 2)
  
  # Skip if no ucos were detected
  skip_if(is.null(cb_res$ucos_details), "No ucos detected in test data")
  
  # Apply both filters
  expect_message(
    result <- get_robust_ucos(cb_res, npc_outcome_cutoff = 0.3, pvalue_cutoff = 0.01),
    "Extracting outcome-specific.*pvalue_cutoff = 0.01.*npc_outcome_cutoff = 0.3"
  )
  
  # Should return a colocboost object
  expect_s3_class(result, "colocboost")
})

test_that("get_robust_ucos handles npc_outcome_cutoff = 0 correctly", {
  
  # Generate test data
  cb_res <- generate_ucos_test_data(output_level = 2)
  
  # Skip if no ucos were detected
  skip_if(is.null(cb_res$ucos_details), "No ucos detected in test data")
  
  # With npc_outcome_cutoff = 0 and no pvalue_cutoff, should return unchanged
  expect_message(
    result <- get_robust_ucos(cb_res, npc_outcome_cutoff = 0),
    "All possible uncolocalized events are reported"
  )
  
  # Should be essentially unchanged
  expect_equal(
    length(result$ucos_details$ucos$ucos_index),
    length(cb_res$ucos_details$ucos$ucos_index)
  )
})

test_that("get_robust_ucos handles missing ucos_details", {
  
  # Generate test data
  cb_res <- generate_ucos_test_data(output_level = 1)  # No ucos_details
  
  # Should warn and return unchanged
  expect_warning(
    result <- get_robust_ucos(cb_res),
    "Please run colocboost model with output_level=2"
  )
  
  # Should return original object
  expect_equal(result, cb_res)
})


test_that("get_robust_ucos validates input object type", {
  
  # Should error with non-colocboost object
  expect_error(
    get_robust_ucos("not_a_colocboost_object"),
    "Input must from colocboost object!"
  )
  
  expect_error(
    get_robust_ucos(list(some = "data")),
    "Input must from colocboost object!"
  )
})

test_that("get_robust_ucos validates pvalue_cutoff range", {
  
  # Generate test data
  cb_res <- generate_ucos_test_data(output_level = 2)
  
  # Skip if no ucos were detected
  skip_if(is.null(cb_res$ucos_details), "No ucos detected in test data")
  
  # pvalue_cutoff > 1 should warn
  expect_warning(
    result <- get_robust_ucos(cb_res, pvalue_cutoff = 1.5),
    "Please check the pvalue cutoff in \\[0,1\\]"
  )
  
  # pvalue_cutoff < 0 should warn
  expect_warning(
    result <- get_robust_ucos(cb_res, pvalue_cutoff = -0.1),
    "Please check the pvalue cutoff in \\[0,1\\]"
  )
})


test_that("get_robust_ucos correctly removes all ucos when all fail cutoff", {
  
  # Generate test data
  cb_res <- generate_ucos_test_data(output_level = 2)
  
  # Skip if no ucos were detected
  skip_if(is.null(cb_res$ucos_details), "No ucos detected in test data")
  
  # Apply impossible cutoff
  result <- get_robust_ucos(cb_res, npc_outcome_cutoff = 1.0)
  
  # Should return colocboost object
  expect_s3_class(result, "colocboost")
  
  # ucos_details should be NULL
  expect_null(result$ucos_details)
})

test_that("get_robust_ucos maintains data structure integrity", {
  
  # Generate test data
  cb_res <- generate_ucos_test_data(output_level = 2)
  
  # Skip if no ucos were detected
  skip_if(is.null(cb_res$ucos_details), "No ucos detected in test data")
  
  # Apply moderate filtering
  result <- get_robust_ucos(cb_res, npc_outcome_cutoff = 0.2)
  
  # Should maintain colocboost structure
  expect_s3_class(result, "colocboost")
  
  # Should have expected fields
  expect_true("data_info" %in% names(result))
  expect_true("model_info" %in% names(result))
  
  # If ucos_details exists, check its structure
  if (!is.null(result$ucos_details)) {
    expect_true("ucos" %in% names(result$ucos_details))
    expect_true("ucos_outcomes" %in% names(result$ucos_details))
    expect_true("ucos_outcomes_npc" %in% names(result$ucos_details))
    expect_true("ucos_weight" %in% names(result$ucos_details))
    expect_true("ucos_purity" %in% names(result$ucos_details))
    
    # Check that all list elements have consistent lengths
    n_ucos <- length(result$ucos_details$ucos$ucos_index)
    expect_equal(length(result$ucos_details$ucos$ucos_variables), n_ucos)
    expect_equal(length(result$ucos_details$ucos_outcomes$outcome_index), n_ucos)
    expect_equal(length(result$ucos_details$ucos_outcomes$outcome_name), n_ucos)
    expect_equal(length(result$ucos_details$ucos_weight), n_ucos)
    expect_equal(nrow(result$ucos_details$ucos_outcomes_npc), n_ucos)
  }
})


test_that("get_robust_ucos preserves names after filtering", {
  
  # Generate test data
  cb_res <- generate_ucos_test_data(output_level = 2)
  
  # Skip if no ucos were detected
  skip_if(is.null(cb_res$ucos_details), "No ucos detected in test data")
  
  # Apply filtering
  result <- get_robust_ucos(cb_res, npc_outcome_cutoff = 0.2)
  
  # If ucos remain, check that names don't contain "remove"
  if (!is.null(result$ucos_details)) {
    ucos_names <- names(result$ucos_details$ucos_outcomes$outcome_index)
    
    # Should not have "remove" in names
    expect_false(any(grepl("^remove", ucos_names)))
    
    # Should have proper "ucos" naming pattern
    expect_true(all(grepl("^ucos", ucos_names)))
  }
})


# Utility function to generate test data with exactly ONE uncolocalized effect
generate_single_ucos_test_data <- function(n = 200, p = 40, seed = 123) {
  set.seed(seed)
  
  # Only 2 traits for simplicity
  L <- 2
  
  # Generate X with LD structure
  sigma <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigma[i, j] <- 0.9^abs(i - j)
    }
  }
  X <- MASS::mvrnorm(n, rep(0, p), sigma)
  colnames(X) <- paste0("SNP", 1:p)
  
  # Generate true effects:
  # - Colocalized effect at SNP10 (affects both traits)
  # - ONE trait-specific effect at SNP25 (affects only trait 2)
  true_beta <- matrix(0, p, L)
  true_beta[10, 1] <- 0.8  # SNP10 affects trait 1
  true_beta[10, 2] <- 0.7  # SNP10 also affects trait 2 (colocalized)
  true_beta[25, 2] <- 0.6  # SNP25 ONLY affects trait 2 (this is the single ucos)
  
  # Generate Y with some noise
  Y <- matrix(0, n, L)
  for (l in 1:L) {
    Y[, l] <- X %*% true_beta[, l] + rnorm(n, 0, 0.8)  # Less noise for clearer signals
  }
  
  # Prepare input for colocboost
  Y_input <- lapply(1:L, function(l) Y[,l])
  X_input <- replicate(L, X, simplify = FALSE)
  
  # Run colocboost with output_level = 2 to get ucos_details
  suppressWarnings({
    result <- colocboost(
      X = X_input, 
      Y = Y_input,
      M = 10,  # More iterations for better detection
      output_level = 2
    )
  })
  
  return(result)
}


test_that("get_robust_ucos handles edge case with single ucos", {
  
  # Generate test data
  cb_res <- generate_single_ucos_test_data()
  
  # Apply lenient filtering (should keep the single ucos)
  result_keep <- get_robust_ucos(cb_res, npc_outcome_cutoff = 0.0)
  
  # Should maintain structure
  expect_s3_class(result_keep, "colocboost")
  expect_true(!is.null(result_keep$ucos_details))
  expect_equal(length(result_keep$ucos_details$ucos$ucos_index), 1)
  
  # Apply strict filtering (should remove the single ucos)
  result_remove <- get_robust_ucos(cb_res, npc_outcome_cutoff = 1.0)
  
  # Should have NULL ucos_details
  expect_null(result_remove$ucos_details)
})


# Helper to generate proper cb_obj structure for testing get_ucos_evidence
# Following test_model.R generate_test_model pattern
generate_test_cb_obj_with_ucos <- function(n = 100, p = 20, L = 2, seed = 42) {
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
  
  # Generate true effects with trait-specific (ucos) effects
  true_beta <- matrix(0, p, L)
  true_beta[5, 1] <- 0.5  # SNP5 affects trait 1
  true_beta[5, 2] <- 0.4  # SNP5 also affects trait 2 (colocalized)
  true_beta[10, 2] <- 0.3 # SNP10 only affects trait 2 (trait-specific)
  
  # Generate Y with some noise
  Y <- matrix(0, n, L)
  for (l in 1:L) {
    Y[, l] <- X %*% true_beta[, l] + rnorm(n, 0, 1)
  }
  
  # Convert Y to list
  Y_list <- list(Y[,1], Y[,2])
  X_list <- list(X, X)
  
  # Run colocboost to get diagnostic_details
  suppressWarnings({
    result <- colocboost(
      X = X_list, 
      Y = Y_list,
      M = 5,
      output_level = 3
    )$diagnostic_details
  })
  
  # Reconstruct cb_data properly following test_model.R pattern
  result$cb_model_para$update_y <- c(1:result$cb_model_para$L)
  Y_list <- lapply(Y_list, as.matrix)
  
  result$cb_data <- colocboost_init_data(
    X = X_list,
    Y = Y_list,
    dict_YX = c(1, 2),
    Z = NULL,
    LD = NULL,
    N_sumstat = NULL,
    dict_sumstatLD = NULL,
    Var_y = NULL,
    SeBhat = NULL,
    keep_variables = lapply(X_list, colnames)
  )
  
  class(result) <- "colocboost"
  result
}


# ============================================================================
# Tests for get_ucos_evidence (internal function)
# ============================================================================

test_that("get_ucos_evidence returns correct structure", {
  
  # Generate proper cb_obj structure following test_model.R pattern
  cb_obj <- generate_test_cb_obj_with_ucos()
  
  # Also run colocboost to get ucos_details
  set.seed(42)
  n <- 100
  p <- 20
  L <- 2
  
  sigma <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigma[i, j] <- 0.9^abs(i - j)
    }
  }
  X <- MASS::mvrnorm(n, rep(0, p), sigma)
  colnames(X) <- paste0("SNP", 1:p)
  
  true_beta <- matrix(0, p, L)
  true_beta[5, 1] <- 0.5
  true_beta[5, 2] <- 0.4
  true_beta[10, 2] <- 0.3
  
  Y <- matrix(0, n, L)
  for (l in 1:L) {
    Y[, l] <- X %*% true_beta[, l] + rnorm(n, 0, 1)
  }
  
  Y_list <- list(Y[,1], Y[,2])
  X_list <- list(X, X)
  
  suppressWarnings({
    cb_res <- colocboost(
      X = X_list, 
      Y = Y_list,
      M = 5,
      output_level = 2
    )
  })
  
  # Skip if no ucos were detected
  skip_if(is.null(cb_res$ucos_details), "No ucos detected in test data")
  
  # Try to access the function from the package namespace
  if (exists("get_ucos_evidence", envir = asNamespace("colocboost"), mode = "function")) {
    get_ucos_evidence <- get("get_ucos_evidence", envir = asNamespace("colocboost"))
    
    # Prepare ucoloc_info from the ucos_details
    ucoloc_info <- list(
      ucos = cb_res$ucos_details$ucos$ucos_index,
      outcome = cb_res$ucos_details$ucos_outcomes$outcome_index,
      outcome_name = cb_res$ucos_details$ucos_outcomes$outcome_name
    )
    
    # Call get_ucos_evidence with the proper cb_obj structure
    result <- get_ucos_evidence(cb_obj, ucoloc_info)
    
    # Check structure
    expect_true(is.data.frame(result))
    
    # Check expected columns
    expected_cols <- c("outcome", "outcomes_index", "relative_logLR", "npc_outcome")
    expect_true(all(expected_cols %in% colnames(result)))
    
    # Check that npc_outcome is in [0, 1]
    expect_true(all(result$npc_outcome >= 0 & result$npc_outcome <= 1))
    
    # Check that relative_logLR is non-negative
    expect_true(all(result$relative_logLR >= 0))
    
  } else {
    skip("get_ucos_evidence not accessible for testing")
  }
})

test_that("get_ucos_evidence handles individual-level data", {
  
  # Generate test data
  set.seed(123)
  n <- 200
  p <- 30
  L <- 2
  
  sigma <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigma[i, j] <- 0.9^abs(i - j)
    }
  }
  X <- MASS::mvrnorm(n, rep(0, p), sigma)
  colnames(X) <- paste0("SNP", 1:p)
  
  true_beta <- matrix(0, p, L)
  true_beta[10, 1] <- 0.7
  true_beta[20, 2] <- 0.6
  
  Y <- matrix(0, n, L)
  for (l in 1:L) {
    Y[, l] <- X %*% true_beta[, l] + rnorm(n, 0, 1)
  }
  
  Y_list <- lapply(1:L, function(l) Y[,l])
  X_list <- replicate(L, X, simplify = FALSE)
  
  # Run colocboost to get ucos_details
  suppressWarnings({
    cb_res <- colocboost(
      X = X_list, 
      Y = Y_list,
      M = 5,
      output_level = 3
    )
  })
  
  skip_if(is.null(cb_res$ucos_details), "No ucos detected in test data")
  
  # Build proper cb_obj structure following test_model.R
  cb_obj <- cb_res$diagnostic_details
  cb_obj$cb_model_para$update_y <- c(1:cb_obj$cb_model_para$L)
  Y_list <- lapply(Y_list, as.matrix)
  
  cb_obj$cb_data <- colocboost_init_data(
    X = X_list,
    Y = Y_list,
    dict_YX = c(1, 2),
    Z = NULL,
    LD = NULL,
    N_sumstat = NULL,
    dict_sumstatLD = NULL,
    Var_y = NULL,
    SeBhat = NULL,
    keep_variables = lapply(X_list, colnames)
  )
  class(cb_obj) <- "colocboost"
  
  # Try to access the function from the package namespace
  if (exists("get_ucos_evidence", envir = asNamespace("colocboost"), mode = "function")) {
    get_ucos_evidence <- get("get_ucos_evidence", envir = asNamespace("colocboost"))
    
    ucoloc_info <- list(
      ucos = cb_res$ucos_details$ucos$ucos_index,
      outcome = cb_res$ucos_details$ucos_outcomes$outcome_index,
      outcome_name = cb_res$ucos_details$ucos_outcomes$outcome_name
    )
    
    # Should work with individual-level data
    expect_error(
      result <- get_ucos_evidence(cb_obj, ucoloc_info),
      NA
    )
    
    expect_true(is.data.frame(result))
    
  } else {
    skip("get_ucos_evidence not accessible for testing")
  }
})

test_that("get_ucos_evidence handles summary statistics data", {
  
  # Generate test data
  set.seed(456)
  n <- 200
  p <- 30
  L <- 2
  
  sigma <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigma[i, j] <- 0.9^abs(i - j)
    }
  }
  X <- MASS::mvrnorm(n, rep(0, p), sigma)
  colnames(X) <- paste0("SNP", 1:p)
  
  true_beta <- matrix(0, p, L)
  true_beta[10, 1] <- 0.7
  true_beta[20, 2] <- 0.6
  
  Y <- matrix(0, n, L)
  for (l in 1:L) {
    Y[, l] <- X %*% true_beta[, l] + rnorm(n, 0, 1)
  }
  
  # Calculate summary statistics
  beta <- matrix(0, p, L)
  se <- matrix(0, p, L)
  for (i in 1:L) {
    for (j in 1:p) {
      fit <- summary(lm(Y[,i] ~ X[,j]))$coef
      if (nrow(fit) == 2) {
        beta[j,i] <- fit[2,1]
        se[j,i] <- fit[2,2]
      }
    }
  }
  
  sumstat_list <- lapply(1:L, function(i) {
    data.frame(
      beta = beta[,i],
      sebeta = se[,i],
      n = n,
      variant = colnames(X)
    )
  })
  
  LD_matrix <- cor(X)
  
  # Run colocboost
  suppressWarnings({
    cb_res <- colocboost(
      sumstat = sumstat_list,
      LD = LD_matrix,
      M = 5,
      output_level = 3
    )
  })
  
  skip_if(is.null(cb_res$ucos_details), "No ucos detected in test data")
  
  # Build proper cb_obj structure
  cb_obj <- cb_res$diagnostic_details
  cb_obj$cb_model_para$update_y <- c(1:cb_obj$cb_model_para$L)
  
  Z_list <- lapply(sumstat_list, function(s) s$beta / s$sebeta)
  
  cb_obj$cb_data <- colocboost_init_data(
    X = NULL,
    Y = NULL,
    dict_YX = NULL,
    Z = Z_list,
    LD = list(LD_matrix),
    N_sumstat = lapply(sumstat_list, function(s) s$n[1]),
    dict_sumstatLD = c(1, 1),
    Var_y = NULL,
    SeBhat = NULL,
    keep_variables = lapply(sumstat_list, function(s) s$variant)
  )
  class(cb_obj) <- "colocboost"
  
  # Try to access the function
  if (exists("get_ucos_evidence", envir = asNamespace("colocboost"), mode = "function")) {
    get_ucos_evidence <- get("get_ucos_evidence", envir = asNamespace("colocboost"))
    
    ucoloc_info <- list(
      ucos = cb_res$ucos_details$ucos$ucos_index,
      outcome = cb_res$ucos_details$ucos_outcomes$outcome_index,
      outcome_name = cb_res$ucos_details$ucos_outcomes$outcome_name
    )
    
    # Should work with summary statistics
    expect_error(
      result <- get_ucos_evidence(cb_obj, ucoloc_info),
      NA
    )
    
    expect_true(is.data.frame(result))
    
  } else {
    skip("get_ucos_evidence not accessible for testing")
  }
})

# ============================================================================
# Integration tests
# ============================================================================

test_that("get_robust_ucos and get_ucos_evidence work together", {
  
  # Generate proper cb_obj
  cb_obj <- generate_test_cb_obj_with_ucos(n = 150, p = 30, L = 2, seed = 555)
  
  # Also generate cb_res for filtering
  set.seed(555)
  n <- 150
  p <- 30
  L <- 2
  
  sigma <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigma[i, j] <- 0.9^abs(i - j)
    }
  }
  X <- MASS::mvrnorm(n, rep(0, p), sigma)
  colnames(X) <- paste0("SNP", 1:p)
  
  true_beta <- matrix(0, p, L)
  true_beta[8, 1] <- 0.6
  true_beta[18, 2] <- 0.5
  
  Y <- matrix(0, n, L)
  for (l in 1:L) {
    Y[, l] <- X %*% true_beta[, l] + rnorm(n, 0, 1)
  }
  
  Y_list <- lapply(1:L, function(l) Y[,l])
  X_list <- replicate(L, X, simplify = FALSE)
  
  suppressWarnings({
    cb_res <- colocboost(
      X = X_list, 
      Y = Y_list,
      M = 5,
      output_level = 2
    )
  })
  
  skip_if(is.null(cb_res$ucos_details), "No ucos detected in test data")
  
  # Apply filtering
  filtered_res <- get_robust_ucos(cb_res, npc_outcome_cutoff = 0.2)
  
  # Check that the result is valid
  expect_s3_class(filtered_res, "colocboost")
  
  # If ucos remain, check that we can extract evidence
  if (!is.null(filtered_res$ucos_details)) {
    
    if (exists("get_ucos_evidence", envir = asNamespace("colocboost"), mode = "function")) {
      get_ucos_evidence <- get("get_ucos_evidence", envir = asNamespace("colocboost"))
      
      ucoloc_info <- list(
        ucos = filtered_res$ucos_details$ucos$ucos_index,
        outcome = filtered_res$ucos_details$ucos_outcomes$outcome_index,
        outcome_name = filtered_res$ucos_details$ucos_outcomes$outcome_name
      )
      
      expect_error(
        evidence <- get_ucos_evidence(cb_obj, ucoloc_info),
        NA
      )
      
      # All npc values should meet the cutoff (or be 0 if removed)
      expect_true(all(evidence$npc_outcome >= 0.2 | evidence$npc_outcome == 0))
      
    } else {
      skip("get_ucos_evidence not accessible for testing")
    }
  }
})

test_that("get_robust_ucos with different cutoffs produces expected ordering", {
  
  # Generate test data
  cb_res <- generate_ucos_test_data(output_level = 2)
  
  # Skip if no ucos
  skip_if(is.null(cb_res$ucos_details), "No ucos detected in test data")
  skip_if(length(cb_res$ucos_details$ucos$ucos_index) < 2, "Need at least 2 ucos for this test")
  
  # Apply progressively stricter cutoffs
  cutoffs <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  n_ucos_remaining <- integer(length(cutoffs))
  
  for (i in seq_along(cutoffs)) {
    result <- get_robust_ucos(cb_res, npc_outcome_cutoff = cutoffs[i])
    n_ucos_remaining[i] <- if(is.null(result$ucos_details)) 0 else length(result$ucos_details$ucos$ucos_index)
  }
  
  # The number of remaining ucos should be monotonically non-increasing
  expect_true(all(diff(n_ucos_remaining) <= 0))
})