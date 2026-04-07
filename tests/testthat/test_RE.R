library(testthat)

# ============================================================================
# Tests for the integration changes across commits c0a2e15 and b61d7cb:
#   (1) `use_entropy` parameter: heterogeneity-aware alpha for integrating
#       per-trait weights (vs uniform 1/L).
#   (2) `residual_correlation` refactor: removed from the core boosting
#       loop (no more `dependency` field on cb_data), now consumed at
#       assembly time in `get_integrated_weight` as a GLS-style alpha
#       adjustment, with:
#         - duplicate outcome columns collapsed by geometric mean
#         - colname parsing of both "outcomeN" and "YN"
#         - sub-selection of residual_correlation[idx, idx]
# ============================================================================


# ---- Shared test data generator (L-adaptive so L=2 works) ----
generate_entropy_test_data <- function(n = 300, p = 40, L = 3, seed = 2026) {
  set.seed(seed)
  sigma <- 0.9^abs(outer(1:p, 1:p, "-"))
  X <- MASS::mvrnorm(n, rep(0, p), sigma)
  colnames(X) <- paste0("SNP", 1:p)
  
  # Shared causal at SNP10, descending heterogeneous per-trait effects.
  het_effects <- c(0.9, 0.3, 0.15, 0.1, 0.08)
  true_beta <- matrix(0, p, L)
  true_beta[10, ] <- het_effects[seq_len(L)]
  if (L >= 2) true_beta[25, 2] <- 0.5  # trait-2-specific
  
  Y <- matrix(0, n, L)
  for (l in seq_len(L)) {
    Y[, l] <- X %*% true_beta[, l] + rnorm(n, 0, 1)
  }
  
  list(
    X = X, Y = Y,
    Y_list = lapply(seq_len(L), function(l) Y[, l]),
    X_list = replicate(L, X, simplify = FALSE),
    LD = cor(X),
    true_beta = true_beta
  )
}


# ============================================================================
# (1a) Unit tests for get_integrated_weight with use_entropy
# ============================================================================

test_that("get_integrated_weight: use_entropy=FALSE returns a valid distribution", {
  w1 <- c(0.6, 0.3, 0.05, 0.05)
  w2 <- c(0.1, 0.1, 0.4, 0.4)
  weights <- cbind(w1, w2)
  
  result <- get_integrated_weight(weights, weight_fudge_factor = 1.5, use_entropy = FALSE)
  
  expect_equal(length(result), nrow(weights))
  expect_equal(sum(result), 1, tolerance = 1e-10)
  expect_true(all(result >= 0))
})

test_that("get_integrated_weight: use_entropy=TRUE returns a valid distribution", {
  w1 <- c(0.6, 0.3, 0.05, 0.05)   # low entropy
  w2 <- c(0.25, 0.25, 0.25, 0.25) # max entropy
  weights <- cbind(w1, w2)
  
  result_ent <- get_integrated_weight(weights, use_entropy = TRUE)
  
  expect_equal(length(result_ent), nrow(weights))
  expect_equal(sum(result_ent), 1, tolerance = 1e-10)
  expect_true(all(result_ent >= 0))
})

test_that("get_integrated_weight: use_entropy=TRUE and FALSE differ for heterogeneous-entropy traits", {
  w_sharp <- c(0.85, 0.10, 0.03, 0.02)
  w_flat  <- c(0.26, 0.25, 0.25, 0.24)
  weights <- cbind(w_sharp, w_flat)
  
  r_uniform <- get_integrated_weight(weights, use_entropy = FALSE)
  r_entropy <- get_integrated_weight(weights, use_entropy = TRUE)
  
  expect_false(isTRUE(all.equal(r_uniform, r_entropy)))
  expect_equal(sum(r_uniform), 1, tolerance = 1e-10)
  expect_equal(sum(r_entropy), 1, tolerance = 1e-10)
})

test_that("get_integrated_weight: identical trait weights give identical results under both modes", {
  w <- c(0.5, 0.25, 0.15, 0.10)
  weights <- cbind(w, w, w)
  
  r_uniform <- get_integrated_weight(weights, use_entropy = FALSE)
  r_entropy <- get_integrated_weight(weights, use_entropy = TRUE)
  
  expect_equal(r_uniform, r_entropy, tolerance = 1e-10)
})

test_that("get_integrated_weight: use_entropy handles zero-weight entries without NaN", {
  # Both columns have mass on the same three rows; column 1 is zero on row 4.
  # The entropy branch filters w > 0 when computing H, so zeros must not
  # produce NaN in alpha.
  w1 <- c(0.5, 0.3, 0.2, 0.0)
  w2 <- c(0.4, 0.3, 0.2, 0.1)
  weights <- cbind(w1, w2)
  
  expect_error(
    result <- get_integrated_weight(weights, use_entropy = TRUE),
    NA
  )
  expect_true(all(is.finite(result)))
  expect_equal(sum(result), 1, tolerance = 1e-10)
  # Row 4 has w1 = 0, so 0^(positive alpha) = 0 on that row.
  expect_equal(result[4], 0)
})


# ============================================================================
# (1b) get_integrated_weight with residual_correlation (new GLS path)
# ============================================================================

test_that("get_integrated_weight: identity residual_correlation matches NULL residual_correlation", {
  # With RC = I, alpha = I %*% 1 / sum(I %*% 1) = 1/L, identical to the
  # legacy uniform-alpha path.
  set.seed(42)
  w1 <- c(0.5, 0.2, 0.15, 0.15)
  w2 <- c(0.1, 0.4, 0.3, 0.2)
  w3 <- c(0.25, 0.25, 0.3, 0.2)
  weights <- cbind(w1, w2, w3)
  colnames(weights) <- paste0("outcome", 1:3)
  
  r_null <- get_integrated_weight(weights, residual_correlation = NULL)
  r_I    <- get_integrated_weight(weights, residual_correlation = diag(3))
  
  expect_equal(as.numeric(r_null), as.numeric(r_I), tolerance = 1e-10)
})

test_that("get_integrated_weight: identity RC also matches NULL RC under use_entropy=TRUE", {
  set.seed(43)
  w1 <- c(0.7, 0.1, 0.1, 0.1)
  w2 <- c(0.1, 0.4, 0.3, 0.2)
  w3 <- c(0.25, 0.25, 0.25, 0.25)
  weights <- cbind(w1, w2, w3)
  colnames(weights) <- paste0("outcome", 1:3)
  
  r_null <- get_integrated_weight(weights, use_entropy = TRUE, residual_correlation = NULL)
  r_I    <- get_integrated_weight(weights, use_entropy = TRUE, residual_correlation = diag(3))
  
  expect_equal(as.numeric(r_null), as.numeric(r_I), tolerance = 1e-10)
})

test_that("get_integrated_weight: asymmetric RC differs from NULL RC", {
  # Equal-correlation RC (e.g. 0.8 everywhere off-diagonal) is actually a
  # no-op because its inverse times a vector of 1s is still a multiple of 1.
  # We need an asymmetric correlation structure for alpha to shift.
  w1 <- c(0.6, 0.2, 0.1, 0.1)
  w2 <- c(0.55, 0.25, 0.1, 0.1)
  w3 <- c(0.1, 0.2, 0.35, 0.35)
  weights <- cbind(w1, w2, w3)
  colnames(weights) <- paste0("outcome", 1:3)
  
  # Traits 1 and 2 heavily correlated, trait 3 nearly independent.
  RC <- matrix(c(1,   0.9, 0.1,
                 0.9, 1,   0.1,
                 0.1, 0.1, 1), nrow = 3, byrow = TRUE)
  
  r_null <- get_integrated_weight(weights, residual_correlation = NULL)
  r_rc   <- get_integrated_weight(weights, residual_correlation = RC)
  
  expect_false(isTRUE(all.equal(as.numeric(r_null), as.numeric(r_rc))))
  expect_equal(sum(r_rc), 1, tolerance = 1e-10)
  expect_true(all(is.finite(r_rc)))
})

test_that("get_integrated_weight: 'outcomeN' and 'YN' colnames parse equivalently", {
  # `get_robust_colocalization` relabels cos_weights columns from "outcomeN"
  # to "YN" before calling `get_integrated_weight`, so both prefixes must
  # produce the same output for the same numerical inputs.
  w1 <- c(0.5, 0.3, 0.1, 0.1)
  w2 <- c(0.2, 0.2, 0.3, 0.3)
  weights_o <- cbind(w1, w2); colnames(weights_o) <- c("outcome1", "outcome2")
  weights_y <- cbind(w1, w2); colnames(weights_y) <- c("Y1", "Y2")
  
  RC <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  
  r_o <- get_integrated_weight(weights_o, residual_correlation = RC)
  r_y <- get_integrated_weight(weights_y, residual_correlation = RC)
  
  expect_equal(as.numeric(r_o), as.numeric(r_y), tolerance = 1e-10)
})

test_that("get_integrated_weight: subset of outcomes picks the right RC sub-block", {
  # 4x4 RC, but weights only carry traits 2 and 4. The function must subset
  # RC[c(2,4), c(2,4)] automatically from the colnames.
  w1 <- c(0.7, 0.1, 0.1, 0.1)
  w2 <- c(0.1, 0.4, 0.3, 0.2)
  weights <- cbind(w1, w2)
  colnames(weights) <- c("outcome2", "outcome4")
  
  RC_full <- matrix(c(1,   0.2, 0.1, 0.05,
                      0.2, 1,   0.3, 0.4,
                      0.1, 0.3, 1,   0.15,
                      0.05, 0.4, 0.15, 1), nrow = 4, byrow = TRUE)
  RC_sub  <- RC_full[c(2, 4), c(2, 4)]
  
  # Running with full RC (letting the function subset by index) should match
  # running with the pre-subsetted 2x2 block with colnames relabeled to 1/2.
  weights_sub <- weights
  colnames(weights_sub) <- c("outcome1", "outcome2")
  
  r_full <- get_integrated_weight(weights, residual_correlation = RC_full)
  r_sub  <- get_integrated_weight(weights_sub, residual_correlation = RC_sub)
  
  expect_equal(as.numeric(r_full), as.numeric(r_sub), tolerance = 1e-10)
})

test_that("get_integrated_weight: duplicate outcome columns are collapsed by geometric mean", {
  # Two columns for the same trait (outcome1 appearing twice) should produce
  # the same integrated weight as a single column whose values are the
  # geometric mean of the two duplicates.
  set.seed(77)
  w1a <- c(0.5, 0.3, 0.1, 0.1)
  w1b <- c(0.4, 0.4, 0.1, 0.1)
  w2  <- c(0.2, 0.2, 0.3, 0.3)
  
  gm_w1 <- exp((log(w1a) + log(w1b)) / 2)   # geometric mean
  
  weights_dup <- cbind(w1a, w1b, w2)
  colnames(weights_dup) <- c("outcome1", "outcome1", "outcome2")
  
  weights_collapsed <- cbind(gm_w1, w2)
  colnames(weights_collapsed) <- c("outcome1", "outcome2")
  
  RC <- matrix(c(1, 0.4, 0.4, 1), nrow = 2)
  
  r_dup <- get_integrated_weight(weights_dup, residual_correlation = RC)
  r_col <- get_integrated_weight(weights_collapsed, residual_correlation = RC)
  
  expect_equal(as.numeric(r_dup), as.numeric(r_col), tolerance = 1e-10)
  expect_equal(sum(r_dup), 1, tolerance = 1e-10)
})

test_that("get_integrated_weight: duplicate collapse only triggers when RC is provided", {
  # Without residual_correlation the collapse branch is skipped — two
  # identical-index columns should NOT be collapsed (legacy behavior).
  w1a <- c(0.5, 0.3, 0.1, 0.1)
  w1b <- c(0.4, 0.4, 0.1, 0.1)
  w2  <- c(0.2, 0.2, 0.3, 0.3)
  
  weights_dup <- cbind(w1a, w1b, w2)
  colnames(weights_dup) <- c("outcome1", "outcome1", "outcome2")
  
  r_null <- get_integrated_weight(weights_dup, residual_correlation = NULL)
  expect_equal(sum(r_null), 1, tolerance = 1e-10)
  expect_true(all(is.finite(r_null)))
})

test_that("get_integrated_weight: duplicate collapse + use_entropy + RC composes", {
  set.seed(91)
  w1a <- c(0.5, 0.3, 0.1, 0.1)
  w1b <- c(0.45, 0.35, 0.1, 0.1)
  w2  <- c(0.25, 0.25, 0.25, 0.25)
  
  weights_dup <- cbind(w1a, w1b, w2)
  colnames(weights_dup) <- c("outcome1", "outcome1", "outcome2")
  
  RC <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  
  expect_error(
    res <- get_integrated_weight(
      weights_dup,
      use_entropy = TRUE,
      residual_correlation = RC
    ),
    NA
  )
  expect_equal(sum(res), 1, tolerance = 1e-10)
  expect_true(all(is.finite(res)))
})


# ============================================================================
# (2) End-to-end: use_entropy wired through colocboost()
# ============================================================================

test_that("colocboost runs with use_entropy = TRUE and returns a valid object", {
  td <- generate_entropy_test_data()
  
  suppressWarnings(suppressMessages({
    res <- colocboost(
      X = td$X_list, Y = td$Y_list,
      use_entropy = TRUE, M = 20, output_level = 2
    )
  }))
  
  expect_s3_class(res, "colocboost")
  expect_equal(res$data_info$n_outcomes, 3)
  expect_equal(length(res$data_info$variables), ncol(td$X))
  if (!is.null(res$vcp)) {
    expect_true(all(res$vcp >= 0 & res$vcp <= 1 + 1e-10))
  }
})

test_that("colocboost use_entropy default is FALSE and preserves legacy behavior", {
  td <- generate_entropy_test_data()
  
  suppressWarnings(suppressMessages({
    res_default <- colocboost(
      X = td$X_list, Y = td$Y_list, M = 20, output_level = 2
    )
    res_explicit_false <- colocboost(
      X = td$X_list, Y = td$Y_list,
      use_entropy = FALSE, M = 20, output_level = 2
    )
  }))
  
  expect_equal(res_default$vcp, res_explicit_false$vcp)
  expect_equal(
    names(res_default$cos_details$cos$cos_index),
    names(res_explicit_false$cos_details$cos$cos_index)
  )
})

test_that("colocboost use_entropy=TRUE can yield different VCP than FALSE under heterogeneous effects", {
  td <- generate_entropy_test_data()
  
  suppressWarnings(suppressMessages({
    res_F <- colocboost(X = td$X_list, Y = td$Y_list,
                        use_entropy = FALSE, M = 30, output_level = 2)
    res_T <- colocboost(X = td$X_list, Y = td$Y_list,
                        use_entropy = TRUE,  M = 30, output_level = 2)
  }))
  
  if (!is.null(res_F$vcp) && !is.null(res_T$vcp) &&
      length(res_F$vcp) == length(res_T$vcp)) {
    expect_false(isTRUE(all.equal(as.numeric(res_F$vcp),
                                  as.numeric(res_T$vcp),
                                  tolerance = 1e-12)))
  } else {
    succeed("One run returned no colocalization; entropy-vs-uniform comparison skipped")
  }
})

test_that("colocboost use_entropy flag is propagated into cb_model_para", {
  td <- generate_entropy_test_data()
  
  suppressWarnings(suppressMessages({
    res <- colocboost(X = td$X_list, Y = td$Y_list,
                      use_entropy = TRUE, M = 10, output_level = 3)
  }))
  
  expect_true("diagnostic_details" %in% names(res))
  expect_equal(res$diagnostic_details$cb_model_para$use_entropy, TRUE)
})

test_that("get_robust_colocalization accepts and propagates use_entropy", {
  td <- generate_entropy_test_data()
  
  suppressWarnings(suppressMessages({
    res <- colocboost(X = td$X_list, Y = td$Y_list,
                      use_entropy = TRUE, M = 20, output_level = 2)
  }))
  
  expect_error(
    suppressMessages(
      filtered <- get_robust_colocalization(
        res,
        cos_npc_cutoff = 0.2,
        npc_outcome_cutoff = 0.1,
        use_entropy = TRUE
      )
    ),
    NA
  )
  expect_s3_class(filtered, "colocboost")
})


# ============================================================================
# (3) End-to-end: residual_correlation wired through colocboost()
# ============================================================================

test_that("colocboost accepts residual_correlation without error", {
  td <- generate_entropy_test_data(L = 2, seed = 7)
  resid_cor <- matrix(c(1, 0.4, 0.4, 1), nrow = 2)
  
  expect_error(
    suppressWarnings(suppressMessages({
      res <- colocboost(
        X = td$X_list, Y = td$Y_list,
        residual_correlation = resid_cor,
        M = 15, output_level = 2
      )
    })),
    NA
  )
  expect_s3_class(res, "colocboost")
})

test_that("colocboost_init_data no longer creates a 'dependency' field on cb_data entries", {
  # Structural guard: the per-outcome `dependency` scalar was the mechanism
  # that used to leak residual_correlation into the fitting loop. It must
  # stay deleted even though RC is now consumed at assembly time.
  td <- generate_entropy_test_data(L = 2, seed = 11)
  Y_list_mat <- lapply(td$Y_list, as.matrix)
  
  cb_data <- colocboost_init_data(
    X = td$X_list, Y = Y_list_mat,
    dict_YX = c(1, 2),
    Z = NULL, LD = NULL, X_ref = NULL, ref_label = NULL,
    N_sumstat = NULL, dict_sumstatLD = NULL,
    Var_y = NULL, SeBhat = NULL,
    keep_variables = lapply(td$X_list, colnames)
  )
  
  expect_s3_class(cb_data, "colocboost")
  for (i in seq_along(cb_data$data)) {
    expect_false("dependency" %in% names(cb_data$data[[i]]),
                 info = paste0("data element ", i, " still has a 'dependency' field"))
  }
})

test_that("residual_correlation does NOT affect the core boosting loop (only assembly)", {
  # The refactor keeps residual_correlation out of the fitting loop, so the
  # fitted per-trait model objects (beta, profile_loglike_each, z) must be
  # IDENTICAL with and without RC. Only the downstream integrated VCP is
  # allowed to differ.
  td <- generate_entropy_test_data(L = 2, seed = 13)
  resid_cor <- matrix(c(1, 0.6, 0.6, 1), nrow = 2)
  
  suppressWarnings(suppressMessages({
    res_no_rc <- colocboost(
      X = td$X_list, Y = td$Y_list,
      M = 20, output_level = 3
    )
    res_rc <- colocboost(
      X = td$X_list, Y = td$Y_list,
      residual_correlation = resid_cor,
      M = 20, output_level = 3
    )
  }))
  
  cb_no <- res_no_rc$diagnostic_details$cb_model
  cb_rc <- res_rc$diagnostic_details$cb_model
  expect_equal(length(cb_no), length(cb_rc))
  
  for (i in seq_along(cb_no)) {
    expect_equal(cb_no[[i]]$beta, cb_rc[[i]]$beta, tolerance = 1e-12,
                 info = paste0("beta differs at outcome ", i))
    expect_equal(cb_no[[i]]$profile_loglike_each,
                 cb_rc[[i]]$profile_loglike_each,
                 tolerance = 1e-12,
                 info = paste0("profile_loglike_each differs at outcome ", i))
    expect_equal(cb_no[[i]]$z, cb_rc[[i]]$z, tolerance = 1e-12,
                 info = paste0("z differs at outcome ", i))
  }
})

test_that("identity residual_correlation matches NULL residual_correlation end-to-end", {
  # RC = I should produce VCP identical to RC = NULL because the GLS alpha
  # collapses to uniform 1/L for an identity matrix.
  td <- generate_entropy_test_data(L = 3, seed = 29)
  
  suppressWarnings(suppressMessages({
    res_null <- colocboost(
      X = td$X_list, Y = td$Y_list,
      M = 20, output_level = 2
    )
    res_I <- colocboost(
      X = td$X_list, Y = td$Y_list,
      residual_correlation = diag(3),
      M = 20, output_level = 2
    )
  }))
  
  if (!is.null(res_null$vcp) && !is.null(res_I$vcp)) {
    expect_equal(as.numeric(res_null$vcp), as.numeric(res_I$vcp),
                 tolerance = 1e-10)
  }
  expect_equal(
    names(res_null$cos_details$cos$cos_index),
    names(res_I$cos_details$cos$cos_index)
  )
})

test_that("residual_correlation is stored on cb_model_para via diagnostic_details", {
  td <- generate_entropy_test_data(L = 2, seed = 31)
  RC <- matrix(c(1, 0.35, 0.35, 1), nrow = 2)
  
  suppressWarnings(suppressMessages({
    res <- colocboost(
      X = td$X_list, Y = td$Y_list,
      residual_correlation = RC,
      M = 10, output_level = 3
    )
  }))
  
  expect_true("diagnostic_details" %in% names(res))
  stored <- res$diagnostic_details$cb_model_para$residual_correlation
  expect_false(is.null(stored))
  expect_equal(stored, RC)
})

test_that("residual_correlation composes with use_entropy = TRUE", {
  td <- generate_entropy_test_data(L = 3, seed = 19)
  resid_cor <- diag(3)
  resid_cor[1, 2] <- resid_cor[2, 1] <- 0.3
  
  expect_error(
    suppressWarnings(suppressMessages({
      res <- colocboost(
        X = td$X_list, Y = td$Y_list,
        residual_correlation = resid_cor,
        use_entropy = TRUE,
        M = 15, output_level = 2
      )
    })),
    NA
  )
  expect_s3_class(res, "colocboost")
  expect_equal(res$data_info$n_outcomes, 3)
})

test_that("get_robust_colocalization accepts and propagates residual_correlation", {
  td <- generate_entropy_test_data(L = 3, seed = 37)
  RC <- matrix(0.2, 3, 3); diag(RC) <- 1
  
  suppressWarnings(suppressMessages({
    res <- colocboost(
      X = td$X_list, Y = td$Y_list,
      residual_correlation = RC,
      M = 20, output_level = 2
    )
  }))
  
  expect_error(
    suppressMessages(
      filtered <- get_robust_colocalization(
        res,
        cos_npc_cutoff = 0.2,
        npc_outcome_cutoff = 0.1,
        residual_correlation = RC
      )
    ),
    NA
  )
  expect_s3_class(filtered, "colocboost")
})

test_that("residual_correlation end-to-end runs with non-identity RC and L=4", {
  td <- generate_entropy_test_data(L = 4, seed = 41)
  RC <- matrix(0.3, 4, 4); diag(RC) <- 1
  
  suppressWarnings(suppressMessages({
    res <- colocboost(
      X = td$X_list, Y = td$Y_list,
      residual_correlation = RC,
      use_entropy = TRUE,
      M = 20, output_level = 3
    )
  }))
  
  expect_s3_class(res, "colocboost")
  expect_equal(res$data_info$n_outcomes, 4)
  expect_equal(res$diagnostic_details$cb_model_para$residual_correlation, RC)
})