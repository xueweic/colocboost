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
    true_beta[5, 1] <- 1  # SNP5 affects the trait
    true_beta[10, 1] <- 1 # SNP10 also affects the trait
  } else {
    # Multi-trait case
    true_beta[5, 1] <- 1  # SNP5 affects trait 1
    true_beta[5, 2] <- 1  # SNP5 also affects trait 2 (colocalized)
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
  expected <- c(1,2)
  expect_equal(result, expected)
  
  # Test with different coverage
  result2 <- w_cs(w, coverage = 0.9)
  expected2 <- c(1,2,3,4)  # First 4 elements cover 90% since 3 and 4 with the same 0.1 weight 
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

library(testthat)

# Test for get_merge_ordered_with_indices function - Core case handling
test_that("get_merge_ordered_with_indices handles conflicting variable orders correctly", {
  
  # Test Case 1: Conflicting orders (cycle detection)
  vec1 <- c("A", "B", "C")
  vec2 <- c("B", "A", "C")  # B comes before A (conflicting with vec1)
  vec3 <- c("C", "A", "B")
  
  vector_list_cycle <- list(vec1, vec2, vec3)
  
  # Should trigger warning about different orders and use priority-based fallback
  expect_warning(
    result_cycle <- get_merge_ordered_with_indices(vector_list_cycle),
    "Variable names in different datasets have different orders"
  )
  
  # Should prioritize the order from earlier datasets
  expect_equal(result_cycle, c("A", "B", "C"))  # Order from vec1 (first dataset)
  
  # Test Case 2: Simple cycle A->B, B->A
  vec_simple1 <- c("A", "B")
  vec_simple2 <- c("B", "A")
  
  expect_warning(
    result_simple <- get_merge_ordered_with_indices(list(vec_simple1, vec_simple2)),
    "Variable names in different datasets have different orders"
  )
  expect_equal(result_simple, c("A", "B"))  # Order from first dataset
  
  # Test Case 3: No conflicts (topological sort should work)
  vec_topo1 <- c("A", "B", "C")
  vec_topo2 <- c("B", "C", "D")  # B->C is consistent
  vec_topo3 <- c("C", "D", "E")  # C->D is consistent
  
  expect_silent(
    result_topo <- get_merge_ordered_with_indices(list(vec_topo1, vec_topo2, vec_topo3))
  )
  expect_equal(result_topo, c("A", "B", "C", "D", "E"))
  
  # Test Case 4: Complex cycle with multiple conflicts
  vec_complex1 <- c("X", "Y", "Z")
  vec_complex2 <- c("Y", "Z", "X")  # Creates Y->Z->X cycle
  vec_complex3 <- c("Z", "X", "Y")  # Continues the cycle
  
  expect_warning(
    result_complex <- get_merge_ordered_with_indices(list(vec_complex1, vec_complex2, vec_complex3)),
    "Variable names in different datasets have different orders"
  )
  expect_equal(result_complex, c("X", "Y", "Z"))  # Order from first dataset
  
  # Test Case 5: Partial conflicts
  vec_partial1 <- c("A", "B", "C", "D")
  vec_partial2 <- c("B", "A", "E", "F")  # B->A conflicts, but E->F is new
  
  expect_warning(
    result_partial <- get_merge_ordered_with_indices(list(vec_partial1, vec_partial2)),
    "Variable names in different datasets have different orders"
  )
  expect_equal(result_partial, c("A", "B", "C", "D", "E", "F"))  # Priority to first dataset
  
  
  # Test Case 7: All same elements but different orders
  vec_same1 <- c("X", "Y", "Z")
  vec_same2 <- c("Z", "Y", "X")
  vec_same3 <- c("Y", "X", "Z")
  
  expect_warning(
    result_same <- get_merge_ordered_with_indices(list(vec_same1, vec_same2, vec_same3)),
    "Variable names in different datasets have different orders"
  )
  expect_equal(result_same, c("X", "Y", "Z"))  # Priority to first dataset
})

# Test for get_merge_ordered_with_indices edge cases
test_that("get_merge_ordered_with_indices handles edge cases", {
  
  # Test with NULL inputs
  expect_equal(get_merge_ordered_with_indices(list()), NULL)
  
  # Test with single vector
  expect_equal(get_merge_ordered_with_indices(list(c("A", "B", "C"))), c("A", "B", "C"))
  
  # Test with identical vectors
  vec <- c("A", "B", "C")
  expect_silent(
    result <- get_merge_ordered_with_indices(list(vec, vec, vec))
  )
  expect_equal(result, c("A", "B", "C"))
  
  # Test with numeric vectors (converted to character)
  num_vec1 <- c(1, 2, 3)
  num_vec2 <- c(2, 3, 4)
  expect_silent(
    result_num <- get_merge_ordered_with_indices(list(num_vec1, num_vec2))
  )
  expect_equal(result_num, c("1", "2", "3", "4"))
})

# Test for missing coverage in colocboost_init_data function
test_that("colocboost_init_data handles complex dictionary mappings", {
  
  # Create test data with mixed individual and summary data
  set.seed(123)
  n <- 100
  p_ind <- 50
  p_sumstat <- 100
  
  # Individual level data
  X_ind <- matrix(rnorm(n * p_ind), n, p_ind)
  colnames(X_ind) <- paste0("SNP_IND_", 1:p_ind)
  Y1 <- matrix(rnorm(n), n, 1)
  Y2 <- matrix(rnorm(n), n, 1)
  
  # Summary statistics
  sumstat1 <- data.frame(
    beta = rnorm(p_sumstat),
    sebeta = abs(rnorm(p_sumstat, 0, 0.1)),
    n = n,
    variant = paste0("SNP_SUM_", 1:p_sumstat)
  )
  
  sumstat2 <- data.frame(
    beta = rnorm(p_sumstat),
    sebeta = abs(rnorm(p_sumstat, 0, 0.1)),
    n = n,
    variant = paste0("SNP_SUM_", 1:p_sumstat)
  )
  
  # LD matrix
  LD_matrix <- diag(p_sumstat)
  colnames(LD_matrix) <- rownames(LD_matrix) <- paste0("SNP_SUM_", 1:p_sumstat)
  
  # Test case: Mixed individual and summary data with complex mappings
  X_list <- list(X_ind)
  Y_list <- list(Y1, Y2)
  Z_list <- list(sumstat1$beta / sumstat1$sebeta, sumstat2$beta / sumstat2$sebeta)
  LD_list <- list(LD_matrix)
  N_sumstat <- list(n, n)
  
  dict_YX <- c(1, 1)  # Both Y outcomes use the same X
  dict_sumstatLD <- c(1, 1)  # Both sumstats use the same LD
  
  keep_variables <- list(
    colnames(X_ind),        # X variables
    sumstat1$variant,       # Sumstat 1 variables  
    sumstat2$variant        # Sumstat 2 variables
  )
  
  # Test focal outcome mapping with mixed data
  if (exists("colocboost_init_data")) {
    # Test focal outcome index that maps to summary statistics
    expect_error({
      cb_data_mixed <- colocboost_init_data(
        X = X_list,
        Y = Y_list,
        dict_YX = dict_YX,
        Z = Z_list,
        LD = LD_list,
        N_sumstat = N_sumstat,
        dict_sumstatLD = dict_sumstatLD,
        Var_y = NULL,
        SeBhat = NULL,
        keep_variables = keep_variables,
        focal_outcome_idx = 3,  # Should map to first summary statistic
        focal_outcome_variables = TRUE,
        overlap_variables = FALSE
      )
      
      # Should use variables from the first summary statistic
      expect_equal(length(cb_data_mixed$variable.names), p_sumstat)
      expect_true(all(grepl("SNP_SUM_", cb_data_mixed$variable.names)))
    }, NA)
  } else {
    skip("colocboost_init_data not directly accessible")
  }
})

# Test for error handling in input validation
test_that("colocboost input validation handles edge cases", {
  
  
  # Test with very few overlapping variables
  sumstat1 <- data.frame(
    beta = rnorm(5),
    sebeta = abs(rnorm(5, 0, 0.1)),
    variant = paste0("SNP_A_", 1:5)
  )
  
  sumstat2 <- data.frame(
    beta = rnorm(5),
    sebeta = abs(rnorm(5, 0, 0.1)),
    variant = paste0("SNP_B_", 1:5)  # No overlap with sumstat1
  )
  
warnings <- testthat::capture_warnings(
  colocboost(sumstat = list(sumstat1, sumstat2))
)
expect_true(any(grepl("No or only 1 overlapping variables were found", warnings)))
  
  # Test with inconsistent effect matrices
  effect_est_bad <- matrix(rnorm(10), 5, 2)
  effect_se_bad <- matrix(abs(rnorm(12)), 6, 2)  # Different dimensions
  
  expect_warning(
    colocboost(effect_est = effect_est_bad, effect_se = effect_se_bad),
    "effect_est and effect_se should be the same dimension!"
  )
})

# Test for LD matrix edge cases
test_that("colocboost handles LD matrix edge cases", {
  
  # Test with singular LD matrix (rank deficient)
  set.seed(123)
  n <- 100
  p <- 10
  
  # Create summary statistics
  sumstat_sing <- data.frame(
    beta = rnorm(p),
    sebeta = abs(rnorm(p, 0, 0.1)),
    n = n,
    variant = paste0("SNP", 1:p)
  )
  
  # Create singular LD matrix (first row = second row)
  LD_singular <- matrix(0.5, p, p)
  diag(LD_singular) <- 1
  LD_singular[2, ] <- LD_singular[1, ]  # Make second row identical to first
  LD_singular[, 2] <- LD_singular[, 1]  # Make second column identical to first
  colnames(LD_singular) <- rownames(LD_singular) <- paste0("SNP", 1:p)
  
  # Should handle singular LD matrix without crashing
  expect_error({
    suppressWarnings(result_singular <- colocboost(
      sumstat = list(sumstat_sing, sumstat_sing),
      LD = LD_singular,
      M = 1
    ))
  }, NA)
  
  # Test with identity LD matrix
  LD_identity <- diag(p)
  colnames(LD_identity) <- rownames(LD_identity) <- paste0("SNP", 1:p)
  
  expect_error({
    suppressWarnings(result_identity <- colocboost(
      sumstat = list(sumstat_sing, sumstat_sing),
      LD = LD_identity,
      M = 1
    ))
  }, NA)
})

# Test for outcome name handling edge cases
test_that("colocboost handles outcome names correctly", {
  
  # Generate basic test data
  set.seed(123)
  X_basic <- matrix(rnorm(100 * 10), 100, 10)
  colnames(X_basic) <- paste0("SNP", 1:10)
  Y_basic <- matrix(rnorm(100 * 3), 100, 3)
  
  Y_basic_list <- list(Y_basic[,1], Y_basic[,2], Y_basic[,3])
  X_basic_list <- list(X_basic, X_basic, X_basic)
  

  warnings <- testthat::capture_warnings(
    result_mismatch <- colocboost(
      X = X_basic_list, 
      Y = Y_basic_list,
      outcome_names = c("Trait1", "Trait2"),  # Only 2 names for 3 traits
      M = 2
    )
  )
  expect_true(any(grepl("There are 3 outcomes inputed, but only 2 provided.*Using default outcome_names as Y1,...,YL", warnings)))
  
  expect_equal(result_mismatch$data_info$outcome_info$outcome_names, c("Y1", "Y2", "Y3"))
  
  # Test with NULL outcome_names
  expect_error({
    suppressWarnings(result_null_names <- colocboost(
      X = X_basic_list, 
      Y = Y_basic_list,
      outcome_names = NULL,
      M = 2
    ))
    expect_equal(result_null_names$data_info$outcome_info$outcome_names, c("Y1", "Y2", "Y3"))
  }, NA)
})


# Test for missing value handling edge cases
test_that("colocboost handles missing values in different scenarios", {
  
  # Generate test data with systematic missing values
  set.seed(123)
  n <- 100
  p <- 15
  X_miss <- matrix(rnorm(n * p), n, p)
  colnames(X_miss) <- paste0("SNP", 1:p)
  Y_miss <- matrix(rnorm(n * 2), n, 2)
  
  # Introduce missing values in specific pattern
  missing_indices <- sample(1:n, 10)
  Y_miss[missing_indices, 1] <- NA
  
  # Different missing indices for second outcome
  missing_indices2 <- sample(1:n, 8)
  Y_miss[missing_indices2, 2] <- NA
  
  Y_miss_list <- list(Y_miss[,1], Y_miss[,2])
  X_miss_list <- list(X_miss, X_miss)
  
  # Should handle missing values by creating separate X matrices
  expect_error({
    suppressWarnings(result_miss <- colocboost(
      X = X_miss_list, 
      Y = Y_miss_list,
      M = 2
    ))
  }, NA)
  
  expect_s3_class(result_miss, "colocboost")
})

# Test for custom parameter combinations
test_that("colocboost handles custom parameter combinations", {
  
  # Generate test data
  set.seed(123)
  X_param <- matrix(rnorm(200 * 25), 200, 25)
  colnames(X_param) <- paste0("SNP", 1:25)
  Y_param <- matrix(rnorm(200 * 2), 200, 2)
  
  Y_param_list <- list(Y_param[,1], Y_param[,2])
  X_param_list <- list(X_param, X_param)
  
  # Test extreme parameter values
  expect_error({
    suppressWarnings(result_extreme <- colocboost(
      X = X_param_list, 
      Y = Y_param_list,
      M = 2,
      lambda = 0.99,           # Very high lambda
      tau = 0.001,             # Very small tau
      learning_rate_init = 0.001,  # Very small learning rate
      jk_equiv_corr = 0.99,    # Very high correlation threshold
      coverage = 0.99,         # Very high coverage
      min_abs_corr = 0.1       # Low correlation threshold
    ))
  }, NA)
  
  expect_s3_class(result_extreme, "colocboost")
})

# Test for plotting edge cases with no colocalization
test_that("colocboost_plot handles results with no colocalization", {
  
  # Create a colocboost result with no colocalization
  set.seed(123)
  X_no_coloc <- matrix(rnorm(100 * 20), 100, 20)
  colnames(X_no_coloc) <- paste0("SNP", 1:20)
  
  # Independent traits (no colocalization)
  Y_no_coloc <- matrix(rnorm(100 * 2) * 0.01, 100, 2)  # Very weak independent signals
  
  Y_no_coloc_list <- list(Y_no_coloc[,1], Y_no_coloc[,2])
  X_no_coloc_list <- list(X_no_coloc, X_no_coloc)
  
  suppressWarnings({
    result_no_coloc <- colocboost(
      X = X_no_coloc_list, 
      Y = Y_no_coloc_list,
      M = 2
    )
  })
  
  # Should plot without error even with no colocalization
  expect_error(suppressWarnings(colocboost_plot(result_no_coloc)), NA)
  
  # Test with different y-axis options
  expect_error(suppressWarnings(colocboost_plot(result_no_coloc, y = "vcp")), NA)
  expect_error(suppressWarnings(colocboost_plot(result_no_coloc, y = "cos_vcp")), NA)
})


# Test for complex coverage scenarios in get_cos
test_that("get_cos handles complex coverage scenarios", {
  
  # Generate test colocboost result
  set.seed(123)
  cb_result_complex <- generate_test_result(n = 200, p = 30, L = 3)
  
  # Test with coverage = 1.0 (should include all variants)
  expect_error({
    result_full_coverage <- get_cos(cb_result_complex, coverage = 1.0)
  }, NA)
  
  # Test with coverage = 0.5 (should include fewer variants)
  expect_error({
    result_half_coverage <- get_cos(cb_result_complex, coverage = 0.5)
  }, NA)
  
  # Generate large genotype matrix for purity testing
  set.seed(456)
  N_large <- 150
  P_large <- 100
  sigma_large <- 0.7^abs(outer(1:P_large, 1:P_large, "-"))
  X_large <- MASS::mvrnorm(N_large, rep(0, P_large), sigma_large)
  colnames(X_large) <- paste0("SNP", 1:P_large)
  
  # Test with very stringent purity thresholds
  expect_error({
    result_stringent <- get_cos(cb_result_complex, 
                               coverage = 0.95, 
                               X = X_large, 
                               min_abs_corr = 0.9,
                               median_abs_corr = 0.95)
  }, NA)
  
  # Should potentially filter out impure CoS
  expect_type(result_stringent, "list")
  expect_named(result_stringent, c("cos", "cos_purity"))
})

# Test for get_hierarchical_clusters with extreme correlation structures
test_that("get_hierarchical_clusters handles extreme correlation structures", {
  
  # Test with perfect block structure
  set.seed(789)
  N <- 100
  P <- 12
  
  # Create perfect block correlation structure
  sigma_block <- matrix(0.1, P, P)
  diag(sigma_block) <- 1
  # Block 1: variables 1-4
  sigma_block[1:4, 1:4] <- 0.95
  # Block 2: variables 5-8  
  sigma_block[5:8, 5:8] <- 0.95
  # Block 3: variables 9-12
  sigma_block[9:12, 9:12] <- 0.95
  
  X_block <- MASS::mvrnorm(N, rep(0, P), sigma_block)
  cormat_block <- cor(X_block)
  
  # Should detect 3 clusters with high min_cluster_corr
  result_block <- get_hierarchical_clusters(cormat_block, min_cluster_corr = 0.9)
  
  expect_equal(ncol(result_block$cluster), 3)
  expect_equal(nrow(result_block$cluster), P)
  
  # Test with all low correlations (should result in many clusters)
  cormat_low <- matrix(0.1, 8, 8)
  diag(cormat_low) <- 1
  
  result_low <- get_hierarchical_clusters(cormat_low, min_cluster_corr = 0.5)
  
  # Should create multiple clusters due to low correlations
  expect_gt(ncol(result_low$cluster), 1)
  
  # Test with mixed correlation structure
  P_mixed <- 6
  sigma_mixed <- matrix(0.2, P_mixed, P_mixed)
  diag(sigma_mixed) <- 1
  # High correlation pair: 1-2
  sigma_mixed[1:2, 1:2] <- 0.9
  # Medium correlation pair: 3-4  
  sigma_mixed[3:4, 3:4] <- 0.6
  # Low correlation singletons: 5, 6
  
  X_mixed <- MASS::mvrnorm(N, rep(0, P_mixed), sigma_mixed)
  cormat_mixed <- cor(X_mixed)
  
  result_mixed <- get_hierarchical_clusters(cormat_mixed, min_cluster_corr = 0.7)
  
  expect_type(result_mixed, "list")
  expect_named(result_mixed, c("cluster", "Q_modularity"))
  expect_equal(nrow(result_mixed$cluster), P_mixed)
})

