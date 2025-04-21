library(testthat)

# Utility function to generate a simple colocboost model object
generate_test_model <- function(n = 100, p = 20, L = 2, seed = 42) {
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
    )$diagnostic_details
  })
  result$cb_model_para$update_y <- c(1:result$cb_model_para$L)
  Y_list <- lapply(Y_list, as.matrix)
  result$cb_data <- colocboost_init_data(
    X = X_list,
    Y = Y_list,
    dict_YX = c(1,2), 
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

# Test for colocboost_init_data
test_that("colocboost_init_data correctly initializes data", {
  
  # Generate test data
  set.seed(42)
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p), n, p)
  colnames(X) <- paste0("SNP", 1:p)
  Y <- matrix(rnorm(n*2), n, 2)
  
  # If the function is exported
  if (exists("colocboost_init_data")) {
    expect_error({
      cb_data <- colocboost_init_data(
        X = list(X),
        Y = list(Y[,1,drop=F], Y[,2,drop=F]),
        dict_YX = c(1,1),
        Z = NULL,
        LD = NULL,
        N_sumstat = NULL,
        dict_sumstatLD = NULL,
        Var_y = NULL,
        SeBhat = NULL,
        keep_variables = list(colnames(X)),
      )
    }, NA)
  } else {
    skip("colocboost_init_data not directly accessible")
  }
})

# Test colocboost_assemble function
test_that("colocboost_assemble processes model results", {
  
  # Generate a test model
  cb_obj <- generate_test_model()
  
  # If the function is exported
  if (exists("colocboost_assemble")) {
    expect_error({
      result <- colocboost_assemble(
        cb_obj,
        coverage = 0.95
      )
    }, NA)
  } else {
    skip("colocboost_assemble not directly accessible")
  }
})

# Test for colocboost_workhorse
test_that("colocboost_workhorse performs boosting iterations", {
  
  # Generate a test model
  cb_obj <- generate_test_model()
  
  # If the workhorse function is exported
  if (exists("colocboost_workhorse")) {
    warnings <- capture_warnings({
      result <- colocboost_workhorse(
        cb_obj$cb_data,
        M = 5  # Small number for testing
      )
    })
    # Check if any of the expected warning patterns are present
    expect_true(
      any(grepl("did not coverage", warnings))
    )
    # Then test the result
    expect_s3_class(result, "colocboost")
  } else {
    skip("colocboost_workhorse not directly accessible")
  }

})
