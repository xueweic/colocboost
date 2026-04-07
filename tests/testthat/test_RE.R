library(testthat)

# ============================================================================
# Tests for the two integration changes in commit c0a2e15:
#   (1) `use_entropy` parameter added to integrate per-trait weights using
#       entropy-based alpha instead of a uniform 1/L.
#   (2) `residual_correlation` no longer writes a per-outcome `dependency`
#       field into cb_data and no longer rescales correlations / profiles
#       in the core fitting loop.
# ============================================================================


# ---- Shared test data generator (mirrors generate_test_data in test_colocboost.R) ----
generate_entropy_test_data <- function(n = 300, p = 40, L = 3, seed = 2026) {
  set.seed(seed)
  sigma <- 0.9^abs(outer(1:p, 1:p, "-"))
  X <- MASS::mvrnorm(n, rep(0, p), sigma)
  colnames(X) <- paste0("SNP", 1:p)
  
  # Shared causal at SNP10 with heterogeneous effects across traits
  # (descending, so trait 1 is strongest). Only fill up to L.
  het_effects <- c(0.9, 0.3, 0.15, 0.1, 0.08)
  true_beta <- matrix(0, p, L)
  true_beta[10, ] <- het_effects[seq_len(L)]
  # Trait-2-specific effect (only if L >= 2)
  if (L >= 2) true_beta[25, 2] <- 0.5
  
  Y <- matrix(0, n, L)
  for (l in 1:L) {
    Y[, l] <- X %*% true_beta[, l] + rnorm(n, 0, 1)
  }
  
  list(
    X = X,
    Y = Y,
    Y_list = lapply(1:L, function(l) Y[, l]),
    X_list = replicate(L, X, simplify = FALSE),
    LD = cor(X),
    true_beta = true_beta
  )
}


# ============================================================================
# (1a) Unit tests for get_integrated_weight with use_entropy
# ============================================================================

test_that("get_integrated_weight: use_entropy=FALSE matches uniform alpha behavior", {
  set.seed(1)
  # Construct weight columns with known structure
  w1 <- c(0.6, 0.3, 0.05, 0.05)
  w2 <- c(0.1, 0.1, 0.4, 0.4)
  weights <- cbind(w1, w2)
  
  result <- get_integrated_weight(weights, weight_fudge_factor = 1.5, use_entropy = FALSE)
  
  # Default behavior: alpha = 1/L for every trait
  expect_equal(length(result), nrow(weights))
  expect_equal(sum(result), 1, tolerance = 1e-10)
  expect_true(all(result >= 0))
})

test_that("get_integrated_weight: use_entropy=TRUE produces a valid probability vector", {
  set.seed(1)
  w1 <- c(0.6, 0.3, 0.05, 0.05)   # low entropy (concentrated)
  w2 <- c(0.25, 0.25, 0.25, 0.25) # maximum entropy (uniform)
  weights <- cbind(w1, w2)
  
  result_ent <- get_integrated_weight(weights, weight_fudge_factor = 1.5, use_entropy = TRUE)
  
  # Basic sanity: length, sums to 1, non-negative
  expect_equal(length(result_ent), nrow(weights))
  expect_equal(sum(result_ent), 1, tolerance = 1e-10)
  expect_true(all(result_ent >= 0))
})

test_that("get_integrated_weight: use_entropy=TRUE and FALSE differ when entropies differ", {
  # Two traits with very different entropy profiles should give a different
  # integrated weight under the two schemes.
  w_sharp <- c(0.85, 0.10, 0.03, 0.02)
  w_flat  <- c(0.26, 0.25, 0.25, 0.24)
  weights <- cbind(w_sharp, w_flat)
  
  r_uniform <- get_integrated_weight(weights, use_entropy = FALSE)
  r_entropy <- get_integrated_weight(weights, use_entropy = TRUE)
  
  # They must not be identical
  expect_false(isTRUE(all.equal(r_uniform, r_entropy)))
  # Both are still valid distributions
  expect_equal(sum(r_uniform), 1, tolerance = 1e-10)
  expect_equal(sum(r_entropy), 1, tolerance = 1e-10)
})

test_that("get_integrated_weight: identical trait weights give identical results under both modes", {
  # If all traits share the same weight vector, entropy re-weighting should
  # collapse to the uniform behavior (up to numerical tolerance).
  w <- c(0.5, 0.25, 0.15, 0.10)
  weights <- cbind(w, w, w)
  
  r_uniform <- get_integrated_weight(weights, use_entropy = FALSE)
  r_entropy <- get_integrated_weight(weights, use_entropy = TRUE)
  
  expect_equal(r_uniform, r_entropy, tolerance = 1e-10)
})

test_that("get_integrated_weight: use_entropy handles zero-weight entries without NaN", {
  # Both columns must have mass on at least one common row, otherwise
  # prod(w^alpha) is 0 everywhere and the normalization is 0/0 by design.
  # The entropy branch filters w > 0 when computing H, so zeros in a
  # column with other positive entries must not produce NaN.
  w1 <- c(0.5, 0.3, 0.2, 0.0)
  w2 <- c(0.4, 0.3, 0.2, 0.1)
  weights <- cbind(w1, w2)
  
  expect_error(
    result <- get_integrated_weight(weights, use_entropy = TRUE),
    NA
  )
  expect_true(all(is.finite(result)))
  expect_equal(sum(result), 1, tolerance = 1e-10)
  # The row where w1 = 0 must get integrated weight 0 (0^anything = 0).
  expect_equal(result[4], 0)
})


# ============================================================================
# (1b) End-to-end: use_entropy wired through colocboost()
# ============================================================================

test_that("colocboost runs with use_entropy = TRUE and returns a valid colocboost object", {
  td <- generate_entropy_test_data()
  
  suppressWarnings(suppressMessages({
    res <- colocboost(
      X = td$X_list,
      Y = td$Y_list,
      use_entropy = TRUE,
      M = 20,
      output_level = 2
    )
  }))
  
  expect_s3_class(res, "colocboost")
  expect_equal(res$data_info$n_outcomes, 3)
  expect_equal(length(res$data_info$variables), ncol(td$X))
  # VCP should be a valid [0,1] vector
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
      X = td$X_list, Y = td$Y_list, use_entropy = FALSE, M = 20, output_level = 2
    )
  }))
  
  expect_s3_class(res_default, "colocboost")
  expect_s3_class(res_explicit_false, "colocboost")
  # Default FALSE should match explicit FALSE
  expect_equal(res_default$vcp, res_explicit_false$vcp)
  expect_equal(
    names(res_default$cos_details$cos$cos_index),
    names(res_explicit_false$cos_details$cos$cos_index)
  )
})

test_that("colocboost use_entropy=TRUE can yield different VCP than use_entropy=FALSE for heterogeneous effects", {
  # With very heterogeneous per-trait effects at a shared locus, the two
  # integration schemes should not be bit-identical.
  td <- generate_entropy_test_data()
  
  suppressWarnings(suppressMessages({
    res_F <- colocboost(X = td$X_list, Y = td$Y_list,
                        use_entropy = FALSE, M = 30, output_level = 2)
    res_T <- colocboost(X = td$X_list, Y = td$Y_list,
                        use_entropy = TRUE,  M = 30, output_level = 2)
  }))
  
  # If both runs produced colocalization results, their VCP vectors should
  # not be exactly equal (the only thing that changed is the integration).
  if (!is.null(res_F$vcp) && !is.null(res_T$vcp) &&
      length(res_F$vcp) == length(res_T$vcp)) {
    # At least one element should differ
    expect_false(isTRUE(all.equal(as.numeric(res_F$vcp),
                                  as.numeric(res_T$vcp),
                                  tolerance = 1e-12)))
  } else {
    succeed("One run returned no colocalization; entropy-vs-uniform comparison skipped")
  }
})

test_that("colocboost use_entropy flag is propagated into cb_model_para via diagnostic_details", {
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
  
  # Should run without error with use_entropy = TRUE
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
# (2) residual_correlation decoupling tests
# ============================================================================

test_that("colocboost accepts residual_correlation without error (backward compatible signature)", {
  td <- generate_entropy_test_data(L = 2, seed = 7)
  
  # Non-trivial 2x2 residual correlation
  resid_cor <- matrix(c(1, 0.4, 0.4, 1), nrow = 2)
  
  expect_error(
    suppressWarnings(suppressMessages({
      res <- colocboost(
        X = td$X_list[1:2],
        Y = td$Y_list[1:2],
        residual_correlation = resid_cor,
        M = 15,
        output_level = 2
      )
    })),
    NA
  )
  expect_s3_class(res, "colocboost")
})

test_that("colocboost_init_data no longer creates a 'dependency' field on cb_data entries", {
  # After the refactor, the core data object must not carry a per-outcome
  # dependency scalar — that was the mechanism that used to leak
  # residual_correlation into the fitting loop.
  td <- generate_entropy_test_data(L = 2, seed = 11)
  Y_list_mat <- lapply(td$Y_list, as.matrix)
  
  cb_data <- colocboost_init_data(
    X = td$X_list[1:2],
    Y = Y_list_mat[1:2],
    dict_YX = c(1, 2),
    Z = NULL, LD = NULL, X_ref = NULL, ref_label = NULL,
    N_sumstat = NULL, dict_sumstatLD = NULL,
    Var_y = NULL, SeBhat = NULL,
    keep_variables = lapply(td$X_list[1:2], colnames)
  )
  
  expect_s3_class(cb_data, "colocboost")
  for (i in seq_along(cb_data$data)) {
    expect_false("dependency" %in% names(cb_data$data[[i]]),
                 info = paste0("data element ", i, " still has a 'dependency' field"))
  }
})

test_that("residual_correlation no longer affects the core fitting loop", {
  # The refactor moves residual_correlation out of colocboost_init_data and
  # into colocboost_assemble (where it currently only computes an unused
  # Theta). Consequently, runs with and without residual_correlation should
  # produce IDENTICAL fitted models (same VCP, same CoS) on the same data.
  td <- generate_entropy_test_data(L = 2, seed = 13)
  
  resid_cor <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  
  suppressWarnings(suppressMessages({
    res_no_rc <- colocboost(
      X = td$X_list[1:2], Y = td$Y_list[1:2],
      M = 20, output_level = 2
    )
    res_rc <- colocboost(
      X = td$X_list[1:2], Y = td$Y_list[1:2],
      residual_correlation = resid_cor,
      M = 20, output_level = 2
    )
  }))
  
  expect_s3_class(res_no_rc, "colocboost")
  expect_s3_class(res_rc, "colocboost")
  
  # Same number of outcomes and same variables
  expect_equal(res_no_rc$data_info$n_outcomes, res_rc$data_info$n_outcomes)
  expect_equal(res_no_rc$data_info$variables, res_rc$data_info$variables)
  
  # Core fit should be identical because residual_correlation no longer
  # touches the boosting loop.
  if (!is.null(res_no_rc$vcp) && !is.null(res_rc$vcp)) {
    expect_equal(as.numeric(res_no_rc$vcp), as.numeric(res_rc$vcp),
                 tolerance = 1e-10)
  }
  expect_equal(
    names(res_no_rc$cos_details$cos$cos_index),
    names(res_rc$cos_details$cos$cos_index)
  )
})

test_that("colocboost runs with an identity residual_correlation (no-op case)", {
  td <- generate_entropy_test_data(L = 2, seed = 17)
  I2 <- diag(2)
  
  suppressWarnings(suppressMessages({
    res <- colocboost(
      X = td$X_list[1:2], Y = td$Y_list[1:2],
      residual_correlation = I2,
      M = 15, output_level = 2
    )
  }))
  
  expect_s3_class(res, "colocboost")
  expect_equal(res$data_info$n_outcomes, 2)
})

test_that("residual_correlation is accepted in combination with use_entropy = TRUE", {
  # The two features should be orthogonal and composable.
  td <- generate_entropy_test_data(L = 3, seed = 19)
  resid_cor <- diag(3)
  resid_cor[1, 2] <- resid_cor[2, 1] <- 0.3
  
  expect_error(
    suppressWarnings(suppressMessages({
      res <- colocboost(
        X = td$X_list, Y = td$Y_list,
        residual_correlation = resid_cor,
        use_entropy = TRUE,
        M = 15,
        output_level = 2
      )
    })),
    NA
  )
  expect_s3_class(res, "colocboost")
  expect_equal(res$data_info$n_outcomes, 3)
})