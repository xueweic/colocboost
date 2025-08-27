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

# Test for dict_keep_variables mapping issue
test_that("colocboost correctly maps focal outcome to keep_variables with dict_keep_variables", {
  
  # Test Case 1: Mixed individual and summary statistics data
  # 5 Y outcomes with 1 superset X + 5 sumstats with 1 superset LD
  set.seed(42)
  n <- 100
  p_x <- 300  # X matrix has 300 variables
  p_sumstat <- 500  # 5th sumstat has 500 variables, others have 300
  
  # Generate superset X matrix
  X_superset <- matrix(rnorm(n * p_x), n, p_x)
  colnames(X_superset) <- paste0("SNP_X_", 1:p_x)
  
  # Generate Y outcomes (all use the same superset X)
  Y_list <- list()
  for (i in 1:5) {
    Y_list[[i]] <- matrix(rnorm(n), n, 1)
  }
  
  # Generate summary statistics
  # First 4 sumstats have 300 variables, 5th has 500 variables
  sumstat_list <- list()
  for (i in 1:4) {
    sumstat_list[[i]] <- data.frame(
      beta = rnorm(p_x),
      sebeta = abs(rnorm(p_x, 0, 0.1)),
      n = n,
      variant = paste0("SNP_S_", 1:p_x)
    )
  }
  
  # 5th sumstat has 500 variables
  sumstat_list[[5]] <- data.frame(
    beta = rnorm(p_sumstat),
    sebeta = abs(rnorm(p_sumstat, 0, 0.1)),
    n = n,
    variant = paste0("SNP_S_", 1:p_sumstat)
  )
  
  # Create superset LD matrix for summary statistics
  LD_superset <- diag(p_sumstat)  # Use diagonal for simplicity
  colnames(LD_superset) <- rownames(LD_superset) <- paste0("SNP_S_", 1:p_sumstat)
  
  # Set up dictionaries
  dict_YX <- c(1, 1, 1, 1, 1)  # All 5 Y outcomes use the same X (index 1)
  dict_sumstatLD <- c(1, 1, 1, 1, 1)  # All 5 sumstats use the same LD (index 1)
  
  # Expected keep_variables structure:
  # Index 1: Variables from X_superset (300 variables)
  # Index 2: Variables from sumstat 1 (300 variables) 
  # Index 3: Variables from sumstat 2 (300 variables)
  # Index 4: Variables from sumstat 3 (300 variables)
  # Index 5: Variables from sumstat 4 (300 variables)
  # Index 6: Variables from sumstat 5 (500 variables)
  keep_variables <- list(
    colnames(X_superset),                    # Index 1: 300 variables
    sumstat_list[[1]]$variant,              # Index 2: 300 variables
    sumstat_list[[2]]$variant,              # Index 3: 300 variables
    sumstat_list[[3]]$variant,              # Index 4: 300 variables
    sumstat_list[[4]]$variant,              # Index 5: 300 variables
    sumstat_list[[5]]$variant               # Index 6: 500 variables
  )
  
  # Test the dictionary creation logic
  if (!is.null(dict_YX) & !is.null(dict_sumstatLD)) {
    dict <- c(dict_YX, max(dict_YX) + dict_sumstatLD)
    n_ind_variable <- max(dict_YX)
    # Create dict_keep_variables: first part maps to X matrices, second part maps to sumstat indices
    dict_keep_variables <- c(dict_YX, 1:length(dict_sumstatLD) + n_ind_variable)
  }
  
  # Expected results:
  # dict = c(1,1,1,1,1,2,2,2,2,2) (maps to data matrices)
  # dict_keep_variables = c(1,1,1,1,1,2,3,4,5,6) (maps to keep_variables indices)
  
  expect_equal(dict, c(1,1,1,1,1,2,2,2,2,2))
  expect_equal(dict_keep_variables, c(1,1,1,1,1,2,3,4,5,6))
  
  # Test focal outcome mapping
  focal_outcome_idx <- 10  # 5th sumstat (the one with 500 variables)
  focal_outcome_variables <- TRUE
  overlap_variables <- FALSE
  
  if (focal_outcome_variables & !is.null(focal_outcome_idx)) {
    if (focal_outcome_idx > length(dict)) {
      stop("Target outcome index is over the total number of outcomes! please check!")
    }
    
    # Use dict_keep_variables to get the correct keep_variables index
    keep_variable_names <- keep_variables[[dict_keep_variables[focal_outcome_idx]]]
    
    # The focal outcome (index 10) should map to keep_variables index 6 (500 variables)
    expect_equal(dict_keep_variables[focal_outcome_idx], 6)
    expect_equal(length(keep_variable_names), 500)
    expect_equal(keep_variable_names, sumstat_list[[5]]$variant)
  }
  
  # Test Case 2: Only summary statistics
  dict_YX_null <- NULL
  dict_sumstatLD_only <- c(1,1,1,1,1)
  
  if (is.null(dict_YX_null) & !is.null(dict_sumstatLD_only)) {
    dict_only_sumstat <- dict_sumstatLD_only
    n_ind_variable_only <- 0
    # Only summary statistics, so dict_keep_variables maps directly to sumstat indices
    dict_keep_variables_only <- 1:length(dict_sumstatLD_only)
  }
  
  expect_equal(dict_keep_variables_only, c(1,2,3,4,5))
  
  # Test focal outcome mapping for sumstat-only case
  focal_outcome_idx_only <- 5  # 5th sumstat
  keep_variable_names_only <- keep_variables[[dict_keep_variables_only[focal_outcome_idx_only] + 1]]  # +1 because no X matrix
  
  expect_equal(dict_keep_variables_only[focal_outcome_idx_only], 5)
  expect_equal(length(keep_variable_names_only), 500)
  
  # Test Case 3: Only individual level data
  dict_YX_only <- c(1,1,1,1,1)
  dict_sumstatLD_null <- NULL
  
  if (!is.null(dict_YX_only) & is.null(dict_sumstatLD_null)) {
    dict_only_ind <- dict_YX_only
    # Only individual level data, so dict_keep_variables is the same as dict_YX
    dict_keep_variables_only_ind <- dict_YX_only
  }
  
  expect_equal(dict_keep_variables_only_ind, c(1,1,1,1,1))
  
  # All outcomes should map to the same X matrix (index 1)
  focal_outcome_idx_ind <- 3
  expect_equal(dict_keep_variables_only_ind[focal_outcome_idx_ind], 1)
  
  # Test Case 4: Test with actual colocboost call (if the function accepts the corrected logic)
  # This would be a integration test to ensure the fix works end-to-end
  
  # Create a minimal working example
  # Skip if colocboost doesn't have the fixed dict_keep_variables logic yet
  if (exists("colocboost_init_data")) {
    expect_error({
      # This should work with the corrected dict_keep_variables mapping
      cb_data_test <- colocboost_init_data(
        X = list(X_superset),
        Y = Y_list,
        dict_YX = dict_YX,
        Z = lapply(sumstat_list, function(s) s$beta / s$sebeta),  # z-scores
        LD = list(LD_superset),
        N_sumstat = lapply(sumstat_list, function(s) s$n[1]),
        dict_sumstatLD = dict_sumstatLD,
        Var_y = NULL,
        SeBhat = NULL,
        keep_variables = keep_variables,
        focal_outcome_idx = 10,  # Should correctly map to 500 variables
        focal_outcome_variables = TRUE,
        overlap_variables = FALSE
      )
      
      # The resulting variable names should be from the 5th sumstat (500 variables)
      expect_equal(length(cb_data_test$variable.names), 500)
    }, NA)
  } else {
    skip("colocboost_init_data not directly accessible for integration test")
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
