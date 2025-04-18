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
  
  # Convert Y to list
  Y_list <- list(Y[,1], Y[,2])
  X_list <- list(X, X)
  
  # Run colocboost with minimal parameters to get a model object
  suppressWarnings({
    result <- colocboost(
      X = X_list, 
      Y = Y_list,
      M = 5,  # Small number of iterations for faster testing
      output_level = 3  # Include full model details
    )
  })
  result
}


# Test colocboost_plot function
test_that("colocboost_plot handles different plot options", {
  
  # Generate a test colocboost results
  cb_res <- generate_test_result()
  
  # Basic plot call
  expect_error(suppressWarnings(colocboost_plot(cb_res)), NA)
  
  # Test with different y-axis values
  expect_error(suppressWarnings(colocboost_plot(cb_res, y = "z_original")), NA)
  
  # Test with different outcome_idx
  expect_error(suppressWarnings(colocboost_plot(cb_res, outcome_idx = 1)), NA)
})

# Test get_cos_summary function
test_that("get_cos_summary handles different parameters", {

  # Generate a test colocboost results
  cb_res <- generate_test_result()
  
  # Basic summary call
  expect_error(get_cos_summary(cb_res), NA)
  
  # With custom outcome names
  expect_error(get_cos_summary(cb_res, outcome_names = c("Trait1", "Trait2")), NA)
  
  # With gene name
  summary_with_gene <- get_cos_summary(cb_res, region_name = "TestGene")
  
  # If summary is not NULL, check for region_name column
  if (!is.null(summary_with_gene)) {
    expect_true("region_name" %in% colnames(summary_with_gene))
    expect_equal(summary_with_gene$region_name[1], "TestGene")
  }
})

# Test for get_strong_colocalization
test_that("get_strong_colocalization filters results correctly", {

  # Generate a test colocboost results
  cb_res <- generate_test_result()
  
  # Basic call
  expect_error(get_strong_colocalization(cb_res), NA)
  
  # With stricter thresholds
  expect_error(get_strong_colocalization(cb_res, cos_npc_cutoff = 0.8), NA)
  
  # With p-value threshold
  expect_error(get_strong_colocalization(cb_res, pvalue_cutoff = 0.05), NA)
})