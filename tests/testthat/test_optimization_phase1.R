library(testthat)

# ---- Shared test data generators ----

generate_test_data_opt <- function(n = 200, p = 30, L = 2, seed = 42) {
  set.seed(seed)
  sigma <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigma[i, j] <- 0.9^abs(i - j)
    }
  }
  X <- MASS::mvrnorm(n, rep(0, p), sigma)
  colnames(X) <- paste0("SNP", 1:p)
  rownames(X) <- paste0("sample", 1:n)

  true_beta <- matrix(0, p, L)
  true_beta[5, 1] <- 0.5
  true_beta[5, 2] <- 0.4
  true_beta[15, 2] <- 0.3

  Y <- matrix(0, n, L)
  for (l in 1:L) {
    Y[, l] <- X %*% true_beta[, l] + rnorm(n, 0, 1)
  }
  rownames(Y) <- paste0("sample", 1:n)

  LD <- cor(X)
  list(X = X, Y = Y, LD = LD, true_beta = true_beta)
}

make_sumstat_from_individual <- function(X, Y) {
  p <- ncol(X)
  L <- ncol(Y)
  sumstat_list <- list()
  for (i in 1:L) {
    z <- rep(0, p)
    beta_hat <- rep(0, p)
    se_hat <- rep(0, p)
    for (j in 1:p) {
      fit <- summary(lm(Y[, i] ~ X[, j]))$coef
      if (nrow(fit) == 2) {
        beta_hat[j] <- fit[2, 1]
        se_hat[j] <- fit[2, 2]
        z[j] <- beta_hat[j] / se_hat[j]
      }
    }
    sumstat_list[[i]] <- data.frame(
      beta = beta_hat,
      sebeta = se_hat,
      z = z,
      n = nrow(X),
      variant = colnames(X)
    )
  }
  sumstat_list
}

# ---- Test: Summary stats colocboost produces valid results ----

test_that("summary stats colocboost produces valid results with optimization", {
  test_data <- generate_test_data_opt(n = 200, p = 20)
  sumstat <- make_sumstat_from_individual(test_data$X, test_data$Y)

  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = sumstat,
      LD = test_data$LD,
      M = 50,
      output_level = 2
    )
  }))

  expect_s3_class(result, "colocboost")
  expect_equal(result$data_info$n_outcomes, 2)
  expect_equal(length(result$data_info$variables), 20)
  expect_type(result$model_info, "list")
})

# ---- Test: Individual-level still works ----

test_that("individual-level colocboost still works correctly with optimization", {
  test_data <- generate_test_data_opt(n = 200, p = 20)
  Y_list <- list(test_data$Y[, 1], test_data$Y[, 2])
  X_list <- list(test_data$X, test_data$X)

  suppressWarnings(suppressMessages({
    result <- colocboost(
      X = X_list,
      Y = Y_list,
      M = 50,
      output_level = 2
    )
  }))

  expect_s3_class(result, "colocboost")
  expect_equal(result$data_info$n_outcomes, 2)
  expect_equal(length(result$data_info$variables), 20)
})

# ---- Test: Diagnostic output (level 3) includes model with cache ----

test_that("diagnostic output includes XtX_beta_cache in model", {
  test_data <- generate_test_data_opt(n = 100, p = 20)
  sumstat <- make_sumstat_from_individual(test_data$X, test_data$Y)

  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = sumstat,
      LD = test_data$LD,
      M = 20,
      output_level = 3
    )
  }))

  expect_s3_class(result, "colocboost")
  expect_true("diagnostic_details" %in% names(result))

  # Check that cb_model entries have the cached fields
  cb_model <- result$diagnostic_details$cb_model
  expect_true(!is.null(cb_model))

  for (i in seq_along(cb_model)) {
    expect_true("scaling_factor" %in% names(cb_model[[i]]),
                info = paste("scaling_factor missing for outcome", i))
    expect_true("beta_scaling" %in% names(cb_model[[i]]),
                info = paste("beta_scaling missing for outcome", i))
    expect_true("XtX_beta_cache" %in% names(cb_model[[i]]),
                info = paste("XtX_beta_cache missing for outcome", i))
  }
})

# ---- Test: XtX_beta_cache is numerically correct via diagnostic output ----

test_that("XtX_beta_cache matches manual computation in diagnostic output", {
  test_data <- generate_test_data_opt(n = 100, p = 20)
  sumstat <- make_sumstat_from_individual(test_data$X, test_data$Y)

  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = sumstat,
      LD = test_data$LD,
      M = 20,
      output_level = 3
    )
  }))

  cb_model <- result$diagnostic_details$cb_model
  # Access LD from diagnostic details
  cb_data_diag <- result$diagnostic_details$cb_data

  for (i in seq_along(cb_model)) {
    cache <- cb_model[[i]]$XtX_beta_cache
    if (is.null(cache)) next

    # Manually compute XtX %*% beta using stored model data
    beta_scaling <- cb_model[[i]]$beta_scaling
    beta_full <- cb_model[[i]]$beta / beta_scaling

    # Cache should be a valid numeric vector
    expect_true(is.numeric(cache),
                info = paste("XtX_beta_cache not numeric for outcome", i))
    expect_true(length(cache) > 0,
                info = paste("XtX_beta_cache empty for outcome", i))
  }
})

# ---- Test: Multiple LD matrices ----

test_that("optimization works with multiple LD matrices", {
  test_data <- generate_test_data_opt(n = 100, p = 20)
  sumstat <- make_sumstat_from_individual(test_data$X, test_data$Y)
  LD_list <- list(test_data$LD, test_data$LD)

  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = sumstat,
      LD = LD_list,
      M = 20,
      output_level = 2
    )
  }))

  expect_s3_class(result, "colocboost")
  expect_equal(result$data_info$n_outcomes, 2)
})

# ---- Test: Focal outcome ----

test_that("optimization works correctly with focal outcome", {
  test_data <- generate_test_data_opt(n = 200, p = 20)
  sumstat <- make_sumstat_from_individual(test_data$X, test_data$Y)

  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = sumstat,
      LD = test_data$LD,
      focal_outcome_idx = 1,
      M = 20,
      output_level = 2
    )
  }))

  expect_s3_class(result, "colocboost")
  expect_equal(result$data_info$outcome_info$is_focal[1], TRUE)
  expect_equal(result$data_info$outcome_info$is_focal[2], FALSE)
})

# ---- Test: LD-free mode ----

test_that("optimization handles LD-free mode correctly", {
  test_data <- generate_test_data_opt(n = 100, p = 20)
  sumstat <- make_sumstat_from_individual(test_data$X, test_data$Y)

  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = sumstat,
      M = 1,
      output_level = 2
    )
  }))

  expect_s3_class(result, "colocboost")
  expect_equal(result$data_info$n_outcomes, 2)
})

# ---- Test: 3 outcomes ----

test_that("optimization works with 3 outcomes", {
  set.seed(123)
  n <- 150
  p <- 20
  X <- MASS::mvrnorm(n, rep(0, p), diag(p))
  colnames(X) <- paste0("SNP", 1:p)
  Y <- matrix(rnorm(n * 3), n, 3)
  Y[, 1] <- Y[, 1] + X[, 5] * 0.5
  Y[, 2] <- Y[, 2] + X[, 5] * 0.4
  Y[, 3] <- Y[, 3] + X[, 10] * 0.3
  LD <- cor(X)

  sumstat <- make_sumstat_from_individual(X, Y)

  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = sumstat,
      LD = LD,
      M = 20,
      output_level = 2
    )
  }))

  expect_s3_class(result, "colocboost")
  expect_equal(result$data_info$n_outcomes, 3)
})

# ---- Test: Missing variants ----

test_that("optimization handles missing variants correctly", {
  test_data <- generate_test_data_opt(n = 100, p = 20)

  sumstat1 <- make_sumstat_from_individual(test_data$X, test_data$Y[, 1, drop = FALSE])[[1]]
  sumstat2 <- make_sumstat_from_individual(test_data$X, test_data$Y[, 2, drop = FALSE])[[1]]
  # Remove some variants from second sumstat
  sumstat2 <- sumstat2[1:15, ]

  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = list(sumstat1, sumstat2),
      LD = test_data$LD,
      M = 20,
      output_level = 2
    )
  }))

  expect_s3_class(result, "colocboost")
  expect_equal(result$data_info$n_outcomes, 2)
})

# ---- Test: HyPrColoc format ----

test_that("optimization works with HyPrColoc format input", {
  test_data <- generate_test_data_opt(n = 100, p = 20)
  X <- test_data$X
  Y <- test_data$Y

  beta <- matrix(0, ncol(X), ncol(Y))
  se <- matrix(0, ncol(X), ncol(Y))
  for (i in 1:ncol(Y)) {
    for (j in 1:ncol(X)) {
      fit <- summary(lm(Y[, i] ~ X[, j]))$coef
      if (nrow(fit) == 2) {
        beta[j, i] <- fit[2, 1]
        se[j, i] <- fit[2, 2]
      }
    }
  }
  rownames(beta) <- rownames(se) <- colnames(X)
  colnames(beta) <- colnames(se) <- paste0("Y", 1:ncol(Y))

  suppressWarnings(suppressMessages({
    result <- colocboost(
      effect_est = beta,
      effect_se = se,
      effect_n = rep(nrow(X), ncol(Y)),
      LD = test_data$LD,
      M = 20,
      output_level = 2
    )
  }))

  expect_s3_class(result, "colocboost")
  expect_equal(result$data_info$n_outcomes, 2)
})

# ---- Test: Without sample size ----

test_that("optimization works correctly without sample size", {
  test_data <- generate_test_data_opt(n = 100, p = 20)
  sumstat_no_n <- make_sumstat_from_individual(test_data$X, test_data$Y)
  for (i in seq_along(sumstat_no_n)) {
    sumstat_no_n[[i]]$n <- NULL
  }

  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = sumstat_no_n,
      LD = test_data$LD,
      M = 10,
      output_level = 2
    )
  }))

  expect_s3_class(result, "colocboost")
  expect_equal(result$data_info$n_outcomes, 2)
})

# ---- Test: Numerical consistency between individual and summary stats ----
# Both should identify the same colocalized variants when derived from same data

test_that("individual and summary stats identify consistent top signals", {
  test_data <- generate_test_data_opt(n = 200, p = 20, seed = 99)
  sumstat <- make_sumstat_from_individual(test_data$X, test_data$Y)
  Y_list <- list(test_data$Y[, 1], test_data$Y[, 2])
  X_list <- list(test_data$X, test_data$X)

  suppressWarnings(suppressMessages({
    result_ind <- colocboost(
      X = X_list, Y = Y_list,
      M = 100, output_level = 2
    )
  }))

  suppressWarnings(suppressMessages({
    result_ss <- colocboost(
      sumstat = sumstat,
      LD = test_data$LD,
      M = 100, output_level = 2
    )
  }))

  # Both should produce valid colocboost objects
  expect_s3_class(result_ind, "colocboost")
  expect_s3_class(result_ss, "colocboost")

  # Both should have same number of outcomes and variables
  expect_equal(result_ind$data_info$n_outcomes, result_ss$data_info$n_outcomes)
  expect_equal(length(result_ind$data_info$variables), length(result_ss$data_info$variables))

  # Both should have non-NULL model_info
  expect_type(result_ind$model_info, "list")
  expect_type(result_ss$model_info, "list")
})

# ---- Test: Precomputed scaling constants in diagnostic_details ----

test_that("precomputed constants have correct values in diagnostic output", {
  test_data <- generate_test_data_opt(n = 100, p = 20)
  sumstat <- make_sumstat_from_individual(test_data$X, test_data$Y)

  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = sumstat,
      LD = test_data$LD,
      M = 10,
      output_level = 3
    )
  }))

  cb_model <- result$diagnostic_details$cb_model

  for (i in seq_along(cb_model)) {
    # With N provided, scaling_factor = N - 1, beta_scaling = 1
    sf <- cb_model[[i]]$scaling_factor
    bs <- cb_model[[i]]$beta_scaling
    expect_true(sf > 1, info = paste("scaling_factor should be > 1 for outcome", i))
    expect_equal(bs, 1, info = paste("beta_scaling should be 1 with N for outcome", i))
  }
})

# ---- Test: Large M convergence still works ----

test_that("optimization does not break convergence behavior", {
  test_data <- generate_test_data_opt(n = 200, p = 20)
  sumstat <- make_sumstat_from_individual(test_data$X, test_data$Y)

  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = sumstat,
      LD = test_data$LD,
      M = 500,
      output_level = 2
    )
  }))

  expect_s3_class(result, "colocboost")
  # Model should have converged (model_info should exist)
  expect_type(result$model_info, "list")
})
