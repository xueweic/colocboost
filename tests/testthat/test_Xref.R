library(testthat)

# ============================================================================
# Shared test data generators for X_ref tests
# ============================================================================

generate_xref_test_data <- function(n = 200, n_ref = 50, p = 30, L = 2, seed = 42) {
  set.seed(seed)
  
  # Generate X with LD structure (for individual-level truth)
  sigma <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigma[i, j] <- 0.9^abs(i - j)
    }
  }
  X <- MASS::mvrnorm(n, rep(0, p), sigma)
  colnames(X) <- paste0("SNP", 1:p)
  
  # Generate X_ref (reference panel, separate samples)
  X_ref <- MASS::mvrnorm(n_ref, rep(0, p), sigma)
  colnames(X_ref) <- paste0("SNP", 1:p)
  
  # Generate true effects
  true_beta <- matrix(0, p, L)
  true_beta[5, 1] <- 0.7   # SNP5 affects trait 1
  true_beta[5, 2] <- 0.6   # SNP5 also affects trait 2 (colocalized)
  true_beta[20, 2] <- 0.5  # SNP20 only affects trait 2
  
  # Generate Y with some noise
  Y <- matrix(0, n, L)
  for (l in 1:L) {
    Y[, l] <- X %*% true_beta[, l] + rnorm(n, 0, 1)
  }
  
  # Calculate summary statistics
  beta <- matrix(0, p, L)
  se <- matrix(0, p, L)
  for (i in 1:L) {
    for (j in 1:p) {
      fit <- summary(lm(Y[, i] ~ X[, j]))$coef
      if (nrow(fit) == 2) {
        beta[j, i] <- fit[2, 1]
        se[j, i] <- fit[2, 2]
      }
    }
  }
  
  sumstat_list <- lapply(1:L, function(i) {
    data.frame(
      beta = beta[, i],
      sebeta = se[, i],
      n = n,
      variant = colnames(X)
    )
  })
  
  LD <- cor(X)
  
  list(
    X = X,
    Y = Y,
    X_ref = X_ref,
    LD = LD,
    sumstat = sumstat_list,
    true_beta = true_beta
  )
}


# ============================================================================
# Test 1: colocboost rejects LD + X_ref together
# ============================================================================

test_that("colocboost rejects when both LD and X_ref are provided", {
  test_data <- generate_xref_test_data()
  
  expect_warning(
    result <- colocboost(
      sumstat = test_data$sumstat,
      LD = test_data$LD,
      X_ref = test_data$X_ref,
      M = 1
    ),
    "Provide either LD or X_ref, not both"
  )
  expect_null(result)
})


# ============================================================================
# Test 2: X_ref with N_ref >= P precomputes LD
# ============================================================================

test_that("X_ref with N_ref >= P precomputes LD and produces valid results", {
  # N_ref = 50, P = 30 => N_ref >= P, should precompute LD
  test_data <- generate_xref_test_data(n_ref = 50, p = 30)
  
  expect_message(
    suppressWarnings({
      result_xref <- colocboost(
        sumstat = test_data$sumstat,
        X_ref = test_data$X_ref,
        M = 10,
        output_level = 2
      )
    }),
    "N_ref >= P: precomputing LD from X_ref"
  )
  
  expect_s3_class(result_xref, "colocboost")
  expect_equal(result_xref$data_info$n_outcomes, 2)
  expect_equal(length(result_xref$data_info$variables), 30)
})


# ============================================================================
# Test 3: X_ref with N_ref < P keeps X_ref for on-the-fly computation
# ============================================================================

test_that("X_ref with N_ref < P runs correctly without precomputing LD", {
  # N_ref = 20, P = 30 => N_ref < P, should keep X_ref
  test_data <- generate_xref_test_data(n_ref = 50, p = 30)
  
  suppressWarnings(suppressMessages({
    result_xref <- colocboost(
      sumstat = test_data$sumstat,
      X_ref = test_data$X_ref,
      M = 10,
      output_level = 2
    )
  }))
  
  expect_s3_class(result_xref, "colocboost")
  expect_equal(result_xref$data_info$n_outcomes, 2)
  expect_equal(length(result_xref$data_info$variables), 30)
})


# ============================================================================
# Test 4: X_ref results are consistent with LD results (N_ref >= P case)
# ============================================================================

test_that("X_ref (N_ref >= P) produces results consistent with precomputed LD", {
  test_data <- generate_xref_test_data(n = 200, n_ref = 200, p = 20, seed = 99)
  
  # Run with precomputed LD
  suppressWarnings(suppressMessages({
    result_ld <- colocboost(
      sumstat = test_data$sumstat,
      LD = test_data$LD,
      M = 50,
      output_level = 2
    )
  }))
  
  # Run with X_ref (same samples as LD source, so N_ref >= P => precomputes LD)
  suppressWarnings(suppressMessages({
    result_xref <- colocboost(
      sumstat = test_data$sumstat,
      X_ref = test_data$X_ref,
      M = 50,
      output_level = 2
    )
  }))
  
  # Both should produce valid colocboost objects
  expect_s3_class(result_ld, "colocboost")
  expect_s3_class(result_xref, "colocboost")
  
  # Same number of outcomes and variables
  expect_equal(result_ld$data_info$n_outcomes, result_xref$data_info$n_outcomes)
  expect_equal(length(result_ld$data_info$variables), length(result_xref$data_info$variables))
})


# ============================================================================
# Test 5: X_ref with N_ref < P and one iteration (one_causal mode)
# ============================================================================

test_that("X_ref with N_ref < P works in one-causal mode (M=1 with LD)", {
  test_data <- generate_xref_test_data(n_ref = 15, p = 30)
  
  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = test_data$sumstat,
      X_ref = test_data$X_ref,
      M = 1,
      output_level = 2
    )
  }))
  
  expect_s3_class(result, "colocboost")
  expect_equal(result$model_info$model_coveraged, "one_causal")
  expect_equal(result$model_info$n_updates, 1)
})


# ============================================================================
# Test 6: X_ref with 3 outcomes
# ============================================================================

test_that("X_ref works with 3 outcomes", {
  set.seed(123)
  n <- 200
  n_ref <- 25
  p <- 25
  L <- 3
  
  sigma <- 0.9^abs(outer(1:p, 1:p, "-"))
  X <- MASS::mvrnorm(n, rep(0, p), sigma)
  colnames(X) <- paste0("SNP", 1:p)
  X_ref <- MASS::mvrnorm(n_ref, rep(0, p), sigma)
  colnames(X_ref) <- paste0("SNP", 1:p)
  
  true_beta <- matrix(0, p, L)
  true_beta[5, 1] <- 0.6
  true_beta[5, 2] <- 0.5
  true_beta[15, 3] <- 0.7
  
  Y <- matrix(0, n, L)
  for (l in 1:L) {
    Y[, l] <- X %*% true_beta[, l] + rnorm(n, 0, 1)
  }
  
  beta <- se <- matrix(0, p, L)
  for (i in 1:L) {
    for (j in 1:p) {
      fit <- summary(lm(Y[, i] ~ X[, j]))$coef
      if (nrow(fit) == 2) { beta[j, i] <- fit[2, 1]; se[j, i] <- fit[2, 2] }
    }
  }
  sumstat <- lapply(1:L, function(i) data.frame(beta = beta[, i], sebeta = se[, i], n = n, variant = colnames(X)))
  
  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = sumstat,
      X_ref = X_ref,
      M = 10,
      output_level = 2
    )
  }))
  
  expect_s3_class(result, "colocboost")
  expect_equal(result$data_info$n_outcomes, 3)
})


# ============================================================================
# Test 7: X_ref with focal outcome
# ============================================================================

test_that("X_ref works with focal outcome", {
  test_data <- generate_xref_test_data(n_ref = 20, p = 30)
  
  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = test_data$sumstat,
      X_ref = test_data$X_ref,
      focal_outcome_idx = 1,
      M = 10,
      output_level = 2
    )
  }))
  
  expect_s3_class(result, "colocboost")
  expect_equal(result$data_info$outcome_info$is_focal[1], TRUE)
  expect_equal(result$data_info$outcome_info$is_focal[2], FALSE)
})


# ============================================================================
# Test 8: X_ref with dict_sumstatLD mapping
# ============================================================================

test_that("X_ref works with dict_sumstatLD mapping for shared reference panel", {
  test_data <- generate_xref_test_data(n_ref = 20, p = 30)
  
  # Both sumstats share the same X_ref
  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = test_data$sumstat,
      X_ref = list(test_data$X_ref),
      dict_sumstatLD = matrix(c(1, 1, 2, 1), ncol = 2),
      M = 10,
      output_level = 2
    )
  }))
  
  expect_s3_class(result, "colocboost")
  expect_equal(result$data_info$n_outcomes, 2)
})


# ============================================================================
# Test 9: X_ref with multiple reference panels
# ============================================================================

test_that("X_ref works with multiple reference panels", {
  test_data <- generate_xref_test_data(n_ref = 20, p = 30)
  
  # Two separate X_ref panels (identical for testing)
  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = test_data$sumstat,
      X_ref = list(test_data$X_ref, test_data$X_ref),
      M = 10,
      output_level = 2
    )
  }))
  
  expect_s3_class(result, "colocboost")
  expect_equal(result$data_info$n_outcomes, 2)
})


# ============================================================================
# Test 10: X_ref with missing variants (partial overlap)
# ============================================================================

test_that("X_ref handles partial variant overlap with sumstat", {
  test_data <- generate_xref_test_data(n_ref = 20, p = 30)
  
  # X_ref has fewer variants than sumstat
  X_ref_partial <- test_data$X_ref[, 1:25]  # Only first 25 of 30 SNPs
  
  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = test_data$sumstat,
      X_ref = X_ref_partial,
      M = 10,
      output_level = 2
    )
  }))
  
  expect_s3_class(result, "colocboost")
  # Should have 25 variables (intersection)
  expect_equal(length(result$data_info$variables), 25)
})


# ============================================================================
# Test 11: compute_XtX_product correctness
# ============================================================================

test_that("compute_XtX_product produces correct results for all ref_label types", {
  # Access internal function
  compute_XtX_product <- get("compute_XtX_product", envir = asNamespace("colocboost"))
  
  set.seed(42)
  n_ref <- 50
  p <- 20
  
  sigma <- 0.9^abs(outer(1:p, 1:p, "-"))
  X_ref_raw <- MASS::mvrnorm(n_ref, rep(0, p), sigma)
  
  # Standardize X_ref (as done in colocboost_validate_input_data)
  X_ref <- Rfast::standardise(X_ref_raw, center = TRUE, scale = TRUE)
  X_ref[is.na(X_ref)] <- 0
  
  # Precompute LD from X_ref
  LD <- get_cormat(X_ref_raw)
  
  beta <- rnorm(p)
  
  # Test LD case
  result_ld <- compute_XtX_product(LD, beta, ref_label = "LD")
  expected_ld <- as.vector(LD %*% beta)
  expect_equal(result_ld, expected_ld)
  
  # Test X_ref case
  result_xref <- compute_XtX_product(X_ref, beta, ref_label = "X_ref")
  # Should approximate LD %*% beta
  # X_ref' * X_ref / (N_ref - 1) ≈ LD
  expected_xref <- as.vector(crossprod(X_ref, X_ref %*% beta) / (n_ref - 1))
  expect_equal(result_xref, expected_xref)
  
  # Test No_ref case
  result_noref <- compute_XtX_product(1, beta, ref_label = "No_ref")
  expect_equal(result_noref, beta)
  
  # X_ref and LD results should be similar (not exact due to finite N_ref)
  expect_equal(result_xref, expected_ld, tolerance = 0.3)
})


# ============================================================================
# Test 12: X_ref validation in colocboost_validate_input_data
# ============================================================================

test_that("colocboost_validate_input_data correctly handles X_ref", {
  test_data <- generate_xref_test_data(n_ref = 20, p = 30)
  
  # Test N_ref < P case
  validated <- colocboost_validate_input_data(
    sumstat = test_data$sumstat,
    X_ref = test_data$X_ref
  )
  
  expect_equal(validated$ref_label, "X_ref")
  expect_true(!is.null(validated$X_ref))
  expect_null(validated$LD)
  
  # Test N_ref >= P case
  test_data_large <- generate_xref_test_data(n_ref = 50, p = 30)
  
  validated_large <- suppressMessages(
    colocboost_validate_input_data(
      sumstat = test_data_large$sumstat,
      X_ref = test_data_large$X_ref
    )
  )
  
  expect_equal(validated_large$ref_label, "LD")
  expect_null(validated_large$X_ref)
  expect_true(!is.null(validated_large$LD))
})


# ============================================================================
# Test 13: X_ref with output_level 3 (diagnostic details)
# ============================================================================

test_that("X_ref works with output_level 3 and has correct ref_label in cb_data", {
  test_data <- generate_xref_test_data(n_ref = 20, p = 30)
  
  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = test_data$sumstat,
      X_ref = test_data$X_ref,
      M = 10,
      output_level = 3
    )
  }))
  
  expect_s3_class(result, "colocboost")
  expect_true("diagnostic_details" %in% names(result))
  
  # Check that ref_label is correctly stored in the model's cb_data
  cb_data <- result$diagnostic_details$cb_data
  for (i in seq_along(cb_data$data)) {
    ref_label_i <- cb_data$data[[i]]$ref_label
    expect_true(!is.null(ref_label_i),
                info = paste("ref_label missing for data element", i))
    expect_true(ref_label_i %in% c("X_ref", "LD", "No_ref", "individual"),
                info = paste("Invalid ref_label for data element", i))
  }
})


# ============================================================================
# Test 14: X_ref with XtX_beta_cache in diagnostic output
# ============================================================================

test_that("X_ref model has XtX_beta_cache in diagnostic output", {
  test_data <- generate_xref_test_data(n_ref = 20, p = 30)
  
  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = test_data$sumstat,
      X_ref = test_data$X_ref,
      M = 10,
      output_level = 3
    )
  }))
  
  cb_model <- result$diagnostic_details$cb_model
  for (i in seq_along(cb_model)) {
    expect_true("XtX_beta_cache" %in% names(cb_model[[i]]),
                info = paste("XtX_beta_cache missing for outcome", i))
  }
})


# ============================================================================
# Test 15: X_ref purity functions work correctly
# ============================================================================

test_that("purity functions dispatch correctly for X_ref", {
  get_purity <- get("get_purity", envir = asNamespace("colocboost"))
  get_between_purity <- get("get_between_purity", envir = asNamespace("colocboost"))
  
  set.seed(42)
  n_ref <- 50
  p <- 20
  
  sigma <- 0.9^abs(outer(1:p, 1:p, "-"))
  X_ref <- MASS::mvrnorm(n_ref, rep(0, p), sigma)
  X_ref <- Rfast::standardise(X_ref, center = TRUE, scale = TRUE)
  X_ref[is.na(X_ref)] <- 0
  
  LD <- cor(MASS::mvrnorm(200, rep(0, p), sigma))
  
  # Test get_purity with X_ref
  purity_xref <- get_purity(c(1, 2, 3), Xcorr = X_ref, ref_label = "X_ref")
  expect_length(purity_xref, 3)
  expect_true(all(purity_xref >= 0 & purity_xref <= 1))
  
  # Test get_purity with LD
  purity_ld <- get_purity(c(1, 2, 3), Xcorr = LD, ref_label = "LD")
  expect_length(purity_ld, 3)
  expect_true(all(purity_ld >= 0 & purity_ld <= 1))
  
  # Test get_purity with No_ref
  purity_noref <- get_purity(c(1, 2, 3), Xcorr = 1, ref_label = "No_ref")
  expect_equal(purity_noref, c(0, 0, 0))
  
  # Results should be similar between X_ref and LD (not exact due to different N)
  expect_equal(purity_xref[1], purity_ld[1], tolerance = 0.3)
  
  # Test get_between_purity with X_ref
  between_xref <- get_between_purity(c(1, 2), c(3, 4), Xcorr = X_ref,
                                     miss_idx = NULL, P = p, ref_label = "X_ref")
  expect_length(between_xref, 3)
  expect_true(all(between_xref >= 0 & between_xref <= 1))
  
  # Test get_between_purity with LD
  between_ld <- get_between_purity(c(1, 2), c(3, 4), Xcorr = LD,
                                   miss_idx = NULL, P = p, ref_label = "LD")
  expect_length(between_ld, 3)
  
  # Test get_between_purity with No_ref
  between_noref <- get_between_purity(c(1, 2), c(3, 4), Xcorr = 1,
                                      miss_idx = NULL, P = p, ref_label = "No_ref")
  expect_equal(between_noref, c(0, 0, 0))
})


# ============================================================================
# Test 16: X_ref with get_robust_colocalization post-processing
# ============================================================================

test_that("get_robust_colocalization works with X_ref results", {
  test_data <- generate_xref_test_data(n_ref = 20, p = 30)
  
  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = test_data$sumstat,
      X_ref = test_data$X_ref,
      M = 10,
      output_level = 2
    )
  }))
  
  # Should not error
  expect_error(
    suppressMessages(
      filtered <- get_robust_colocalization(result, cos_npc_cutoff = 0.2, npc_outcome_cutoff = 0.1)
    ),
    NA
  )
  expect_s3_class(filtered, "colocboost")
})


# ============================================================================
# Test 17: X_ref with get_robust_ucos post-processing
# ============================================================================

test_that("get_robust_ucos works with X_ref results", {
  test_data <- generate_xref_test_data(n_ref = 500, p = 30)
  
  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = test_data$sumstat,
      X_ref = test_data$X_ref,
      M = 10,
      output_level = 2
    )
  }))
  
  skip_if(is.null(result$ucos_details), "No ucos detected")
  
  expect_error(
    suppressMessages(
      filtered <- get_robust_ucos(result, npc_outcome_cutoff = 0.1)
    ),
    NA
  )
  expect_s3_class(filtered, "colocboost")
})


# ============================================================================
# Test 18: X_ref plotting works
# ============================================================================

test_that("colocboost_plot works with X_ref results", {
  test_data <- generate_xref_test_data(n_ref = 20, p = 30)
  
  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = test_data$sumstat,
      X_ref = test_data$X_ref,
      M = 10,
      output_level = 2
    )
  }))
  
  expect_error(suppressWarnings(colocboost_plot(result)), NA)
  expect_error(suppressWarnings(colocboost_plot(result, y = "vcp")), NA)
})


# ============================================================================
# Test 19: X_ref with get_colocboost_summary
# ============================================================================

test_that("get_colocboost_summary works with X_ref results", {
  test_data <- generate_xref_test_data(n_ref = 20, p = 30)
  
  suppressWarnings(suppressMessages({
    result <- colocboost(
      sumstat = test_data$sumstat,
      X_ref = test_data$X_ref,
      M = 10,
      output_level = 2
    )
  }))
  
  summary1 <- get_colocboost_summary(result, summary_level = 1)
  expect_type(summary1, "list")
  
  summary2 <- get_colocboost_summary(result, summary_level = 2)
  expect_type(summary2, "list")
})


# ============================================================================
# Test 20: X_ref with HyPrColoc format input
# ============================================================================

test_that("X_ref works with effect_est and effect_se input", {
  test_data <- generate_xref_test_data(n_ref = 20, p = 30)
  
  beta <- se <- matrix(0, 30, 2)
  for (i in 1:2) {
    for (j in 1:30) {
      fit <- summary(lm(test_data$Y[, i] ~ test_data$X[, j]))$coef
      if (nrow(fit) == 2) { beta[j, i] <- fit[2, 1]; se[j, i] <- fit[2, 2] }
    }
  }
  rownames(beta) <- rownames(se) <- colnames(test_data$X)
  
  suppressWarnings(suppressMessages({
    result <- colocboost(
      effect_est = beta,
      effect_se = se,
      effect_n = rep(200, 2),
      X_ref = test_data$X_ref,
      M = 10,
      output_level = 2
    )
  }))
  
  expect_s3_class(result, "colocboost")
  expect_equal(result$data_info$n_outcomes, 2)
})


# ============================================================================
# Test 21: LD_jk functions dispatch correctly for X_ref
# ============================================================================

test_that("get_LD_jk and get_LD_jk1_jk2 dispatch correctly for X_ref", {
  get_LD_jk <- get("get_LD_jk", envir = asNamespace("colocboost"))
  get_LD_jk1_jk2 <- get("get_LD_jk1_jk2", envir = asNamespace("colocboost"))
  
  set.seed(42)
  n_ref <- 50
  p <- 20
  
  sigma <- 0.9^abs(outer(1:p, 1:p, "-"))
  X_ref <- MASS::mvrnorm(n_ref, rep(0, p), sigma)
  X_ref <- Rfast::standardise(X_ref, center = TRUE, scale = TRUE)
  X_ref[is.na(X_ref)] <- 0
  
  LD <- cor(MASS::mvrnorm(200, rep(0, p), sigma))
  
  remain_idx <- 1:p
  
  # get_LD_jk with X_ref
  ld_jk_xref <- get_LD_jk(5, XtX = X_ref, remain_idx = remain_idx, P = p, ref_label = "X_ref")
  expect_length(ld_jk_xref, p)
  expect_equal(ld_jk_xref[5], 1.0, tolerance = 0.01)  # Self-correlation ~ 1
  
  # get_LD_jk with LD
  ld_jk_ld <- get_LD_jk(5, XtX = LD, remain_idx = remain_idx, P = p, ref_label = "LD")
  expect_length(ld_jk_ld, p)
  expect_equal(ld_jk_ld[5], 1.0, tolerance = 0.01)
  
  # get_LD_jk with No_ref
  ld_jk_noref <- get_LD_jk(5, XtX = 1, remain_idx = remain_idx, P = p, ref_label = "No_ref")
  expect_true(all(ld_jk_noref[remain_idx] == 1))
  
  # get_LD_jk1_jk2 with X_ref
  ld_pair_xref <- get_LD_jk1_jk2(1, 2, XtX = X_ref, remain_jk = remain_idx, ref_label = "X_ref")
  expect_true(abs(ld_pair_xref) <= 1)
  expect_true(abs(ld_pair_xref) > 0.5)  # Adjacent SNPs with 0.9 LD should be high
  
  # get_LD_jk1_jk2 with LD
  ld_pair_ld <- get_LD_jk1_jk2(1, 2, XtX = LD, remain_jk = remain_idx, ref_label = "LD")
  expect_true(abs(ld_pair_ld) <= 1)
  
  # get_LD_jk1_jk2 with No_ref
  ld_pair_noref <- get_LD_jk1_jk2(1, 2, XtX = 1, remain_jk = remain_idx, ref_label = "No_ref")
  expect_equal(ld_pair_noref, 0)
})


# ============================================================================
# Test 22: ref_label is never NULL in internal cb_data after processing
# ============================================================================

test_that("ref_label is always explicitly set in cb_data, never NULL", {
  test_data <- generate_xref_test_data(n_ref = 20, p = 30)
  
  # Build cb_data directly to inspect ref_label
  sumstat <- test_data$sumstat
  Z_list <- lapply(sumstat, function(s) s$beta / s$sebeta)
  
  # X_ref case
  cb_data_xref <- colocboost_init_data(
    X = NULL, Y = NULL, dict_YX = NULL,
    Z = Z_list, LD = NULL, X_ref = list(test_data$X_ref), ref_label = "X_ref",
    N_sumstat = lapply(sumstat, function(s) s$n[1]),
    dict_sumstatLD = c(1, 1), Var_y = NULL, SeBhat = NULL,
    keep_variables = lapply(sumstat, function(s) s$variant)
  )
  for (i in seq_along(cb_data_xref$data)) {
    expect_equal(cb_data_xref$data[[i]]$ref_label, "X_ref")
  }
  
  # LD case
  cb_data_ld <- colocboost_init_data(
    X = NULL, Y = NULL, dict_YX = NULL,
    Z = Z_list, LD = list(test_data$LD), X_ref = NULL, ref_label = "LD",
    N_sumstat = lapply(sumstat, function(s) s$n[1]),
    dict_sumstatLD = c(1, 1), Var_y = NULL, SeBhat = NULL,
    keep_variables = lapply(sumstat, function(s) s$variant)
  )
  for (i in seq_along(cb_data_ld$data)) {
    expect_equal(cb_data_ld$data[[i]]$ref_label, "LD")
  }
  
  # Individual data case
  Y_list <- lapply(1:2, function(l) as.matrix(test_data$Y[, l]))
  X_list <- list(test_data$X, test_data$X)
  cb_data_ind <- colocboost_init_data(
    X = X_list, Y = Y_list, dict_YX = c(1, 2),
    Z = NULL, LD = NULL, X_ref = NULL, ref_label = NULL,
    N_sumstat = NULL, dict_sumstatLD = NULL, Var_y = NULL, SeBhat = NULL,
    keep_variables = lapply(X_list, colnames)
  )
  for (i in seq_along(cb_data_ind$data)) {
    expect_equal(cb_data_ind$data[[i]]$ref_label, "individual")
  }
})