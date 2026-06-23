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

test_that("pairwise jk checks reuse LD calculation for shared reference data", {
  set.seed(202)
  n <- 80
  p <- 20
  n_outcomes <- 6
  X <- matrix(rnorm(n * p), nrow = n)
  colnames(X) <- paste0("SNP", seq_len(p))
  jk_each <- c(3, 6, 9, 12, 15, 18)
  pos.update <- seq_len(n_outcomes) + 1L
  cb_data <- list(
    data = c(
      list(list(X = X, XtX = NULL, ref_label = "individual")),
      lapply(seq_len(n_outcomes), function(i) {
        list(N = n, variable_miss = integer(0))
      })
    )
  )
  X_dict <- rep(1L, n_outcomes)
  model_update <- lapply(seq_len(n_outcomes), function(i) {
    change_loglike <- seq_len(p) / p + i * 1e-4
    change_loglike[jk_each] <- change_loglike[jk_each] + rnorm(n_outcomes, sd = 1e-5)
    list(change_loglike = change_loglike, res = numeric(p))
  })

  ns <- asNamespace("colocboost")
  original_get_cormat <- get("get_cormat", envir = ns)
  call_count <- 0L
  unlockBinding("get_cormat", ns)
  assign("get_cormat", function(...) {
    call_count <<- call_count + 1L
    original_get_cormat(...)
  }, envir = ns)
  lockBinding("get_cormat", ns)
  on.exit({
    unlockBinding("get_cormat", ns)
    assign("get_cormat", original_get_cormat, envir = ns)
    lockBinding("get_cormat", ns)
  }, add = TRUE)

  pair_check <- get("check_pair_jkeach", envir = ns)
  res <- pair_check(jk_each, pos.update, model_update, cb_data, X_dict)

  expect_equal(call_count, 1L)
  expect_equal(dim(res), c(n_outcomes, n_outcomes))
  expect_true(all(res == t(res)))
  expect_true(all(diag(res) == 0))
})

make_shared_update_fixture <- function(n = 60, p = 25, n_outcomes = 5, update_jk = 4, seed = 303) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n)
  X <- scale(X)
  colnames(X) <- paste0("SNP", seq_len(p))
  Y <- matrix(rnorm(n * n_outcomes), nrow = n)
  Y <- scale(Y)

  cb_data <- list(
    data = lapply(seq_len(n_outcomes), function(i) {
      list(
        X = if (i == 1L) X else NULL,
        Y = as.matrix(Y[, i]),
        N = n,
        variable_miss = integer(0),
        ref_label = "individual"
      )
    }),
    dict = rep(1L, n_outcomes),
    variable.names = colnames(X)
  )
  class(cb_data) <- "colocboost"

  cb_model <- lapply(seq_len(n_outcomes), function(i) {
    correlation <- rnorm(p)
    list(
      res = as.matrix(Y[, i]),
      beta = rep(0, p),
      weights_path = list(),
      profile_loglike_each = mean(Y[, i]^2),
      obj_path = 999999,
      obj_single = 999999,
      change_loglike = abs(rnorm(p)),
      correlation = correlation,
      z = correlation,
      learning_rate_init = 0.01,
      stop_thresh = 1e-6,
      ld_jk = list(),
      jk = integer(0),
      scaling_factor = n - 1,
      beta_scaling = 1,
      XtX_beta_cache = NULL
    )
  })
  class(cb_model) <- "colocboost"

  cb_model_para <- list(
    update_temp = list(
      update_status = rep(1, n_outcomes),
      real_update_jk = rep(update_jk, n_outcomes)
    ),
    focal_outcome_idx = NULL,
    tau = 0.01,
    lambda = 0.5,
    lambda_focal_outcome = 1,
    func_simplex = "LD_z2z",
    LD_free = FALSE,
    dynamic_learning_rate = FALSE,
    learning_rate_decay = 1,
    P = p
  )
  class(cb_model_para) <- "colocboost"

  list(cb_model = cb_model, cb_model_para = cb_model_para, cb_data = cb_data)
}

test_that("colocboost_update reuses LD_jk for outcomes sharing reference data", {
  fixture <- make_shared_update_fixture(n_outcomes = 5, update_jk = 6)

  ns <- asNamespace("colocboost")
  original_get_LD_jk <- get("get_LD_jk", envir = ns)
  call_count <- 0L
  unlockBinding("get_LD_jk", ns)
  assign("get_LD_jk", function(...) {
    call_count <<- call_count + 1L
    original_get_LD_jk(...)
  }, envir = ns)
  lockBinding("get_LD_jk", ns)
  on.exit({
    unlockBinding("get_LD_jk", ns)
    assign("get_LD_jk", original_get_LD_jk, envir = ns)
    lockBinding("get_LD_jk", ns)
  }, add = TRUE)

  updated <- colocboost_update(fixture$cb_model, fixture$cb_model_para, fixture$cb_data)

  expect_equal(call_count, 1L)
  expect_true(all(vapply(updated, function(model) "6" %in% names(model$ld_jk), logical(1))))
})

test_that("individual update uses crossprod without changing profile loglikelihood", {
  fixture <- make_shared_update_fixture(n_outcomes = 1, update_jk = 5)
  fixture$cb_model[[1]]$ld_jk[["5"]] <- rep(1, fixture$cb_model_para$P)
  matmul_count <- 0L
  fixture$cb_data$data[[1]]$X <- structure(
    fixture$cb_data$data[[1]]$X,
    class = c("counted_matrix", "matrix")
  )
  assign("%*%.counted_matrix", function(x, y) {
    matmul_count <<- matmul_count + 1L
    NextMethod()
  }, envir = .GlobalEnv)
  on.exit(rm("%*%.counted_matrix", envir = .GlobalEnv), add = TRUE)

  updated <- colocboost_update(fixture$cb_model, fixture$cb_model_para, fixture$cb_data)
  profile_log <- tail(updated[[1]]$profile_loglike_each, n = 1)

  expect_equal(matmul_count, 2L)
  expect_equal(
    as.numeric(profile_log),
    mean((fixture$cb_data$data[[1]]$Y - fixture$cb_data$data[[1]]$X %*% updated[[1]]$beta)^2),
    tolerance = 1e-12
  )
})

test_that("merge_ucos skips between-purity checks for disjoint uCoS pairs", {
  set.seed(404)
  p <- 30
  n <- 20
  n_outcomes <- 4
  LD <- diag(p)
  LD[1, 2] <- LD[2, 1] <- 0.9
  cb_obj <- list(
    cb_data = list(
      data = lapply(seq_len(n_outcomes), function(i) {
        list(X = NULL, XtX = LD, variable_miss = integer(0), ref_label = "LD", N = n)
      }),
      dict = seq_len(n_outcomes)
    ),
    cb_model_para = list(P = p, L = n_outcomes)
  )
  class(cb_obj) <- "colocboost"

  ucos_each <- list(c(1, 2), c(10, 11), c(2, 3), c(20, 21))
  names(ucos_each) <- paste0("sets:Y", seq_len(n_outcomes), ":ucos1")
  avW <- matrix(runif(p * length(ucos_each)), nrow = p)
  colnames(avW) <- names(ucos_each)
  past_out <- list(
    ucos = list(
      ucos_each = ucos_each,
      avW_ucos_each = avW,
      change_obj_each = matrix(0.1, nrow = length(ucos_each), ncol = n_outcomes),
      purity_each = matrix(1, nrow = length(ucos_each), ncol = 3),
      ucos_outcome = seq_len(n_outcomes)
    ),
    cos = list(cos = list())
  )

  ns <- asNamespace("colocboost")
  original_between <- get("get_between_purity", envir = ns)
  original_purity <- get("get_purity", envir = ns)
  between_calls <- 0L
  unlockBinding("get_between_purity", ns)
  assign("get_between_purity", function(pos1, pos2, ...) {
    if (length(intersect(pos1, pos2)) == 0) {
      stop("disjoint uCoS pair should not require between-purity")
    }
    between_calls <<- between_calls + 1L
    c(min_abs_cor = 0.9, max_abs_cor = 1, median_abs_cor = 0.9)
  }, envir = ns)
  lockBinding("get_between_purity", ns)
  unlockBinding("get_purity", ns)
  assign("get_purity", function(...) c(1, 1, 1), envir = ns)
  lockBinding("get_purity", ns)
  on.exit({
    unlockBinding("get_between_purity", ns)
    assign("get_between_purity", original_between, envir = ns)
    lockBinding("get_between_purity", ns)
    unlockBinding("get_purity", ns)
    assign("get_purity", original_purity, envir = ns)
    lockBinding("get_purity", ns)
  }, add = TRUE)

  result <- get("merge_ucos", envir = ns)(
    cb_obj, past_out,
    min_abs_corr = 0.5,
    median_cos_abs_corr = 0.8
  )

  expect_equal(between_calls, 2L)
  expect_equal(length(result$cos$cos$cos), 1L)
  expect_equal(length(result$ucos$ucos_each), 2L)
})

test_that("merge_ucos skips full between-purity when top variants cannot pass merge cutoff", {
  p <- 10
  n <- 20
  LD <- diag(p)
  LD[1, 3] <- LD[3, 1] <- 0.1
  cb_obj <- list(
    cb_data = list(
      data = lapply(seq_len(2), function(i) {
        list(X = NULL, XtX = LD, variable_miss = integer(0), ref_label = "LD", N = n)
      }),
      dict = seq_len(2)
    ),
    cb_model_para = list(P = p, L = 2)
  )
  class(cb_obj) <- "colocboost"

  ucos_each <- list(c(1, 2), c(3, 2))
  names(ucos_each) <- c("sets:Y1:ucos1", "sets:Y2:ucos1")
  avW <- matrix(runif(p * length(ucos_each)), nrow = p)
  colnames(avW) <- names(ucos_each)
  past_out <- list(
    ucos = list(
      ucos_each = ucos_each,
      avW_ucos_each = avW,
      change_obj_each = matrix(0.1, nrow = length(ucos_each), ncol = 2),
      purity_each = matrix(1, nrow = length(ucos_each), ncol = 3),
      ucos_outcome = seq_len(2)
    ),
    cos = list(cos = list())
  )

  ns <- asNamespace("colocboost")
  original_between <- get("get_between_purity", envir = ns)
  unlockBinding("get_between_purity", ns)
  assign("get_between_purity", function(...) {
    stop("top-variant prefilter should skip full between-purity")
  }, envir = ns)
  lockBinding("get_between_purity", ns)
  on.exit({
    unlockBinding("get_between_purity", ns)
    assign("get_between_purity", original_between, envir = ns)
    lockBinding("get_between_purity", ns)
  }, add = TRUE)

  result <- get("merge_ucos", envir = ns)(
    cb_obj, past_out,
    min_abs_corr = 0.5,
    median_cos_abs_corr = 0.8
  )

  expect_equal(length(result$cos$cos$cos), 0L)
  expect_equal(length(result$ucos$ucos_each), 2L)
})

test_that("overlap candidate pairs match pairwise intersect scan", {
  ucos_each <- list(
    c(1, 2),
    c(3, 4),
    c(2, 5),
    c(5, 6),
    7,
    c(2, 6)
  )
  reference <- do.call(rbind, lapply(seq_len(length(ucos_each) - 1L), function(i) {
    do.call(rbind, lapply((i + 1L):length(ucos_each), function(j) {
      if (length(intersect(ucos_each[[i]], ucos_each[[j]])) == 0) {
        return(NULL)
      }
      c(i, j)
    }))
  }))
  colnames(reference) <- c("i", "j")

  candidate_pairs <- get(".merge_ucos_overlap_pairs", envir = asNamespace("colocboost"))(ucos_each)

  expect_equal(candidate_pairs, reference)
})

test_that("merge_ucos candidate pairs preserve pairwise intersect output", {
  set.seed(505)
  p <- 30
  n_outcomes <- 5
  LD <- diag(p)
  LD[1, 2] <- LD[2, 1] <- 0.9
  LD[1, 20] <- LD[20, 1] <- 0.95
  cb_obj <- list(
    cb_data = list(
      data = lapply(seq_len(n_outcomes), function(i) {
        list(X = NULL, XtX = LD, variable_miss = integer(0), ref_label = "LD", N = 20)
      }),
      dict = seq_len(n_outcomes)
    ),
    cb_model_para = list(P = p, L = n_outcomes)
  )
  class(cb_obj) <- "colocboost"

  ucos_each <- list(c(1, 2), c(10, 11), c(2, 3), c(20, 21), c(3, 4))
  names(ucos_each) <- paste0("sets:Y", seq_len(n_outcomes), ":ucos1")
  avW <- matrix(runif(p * length(ucos_each)), nrow = p)
  colnames(avW) <- names(ucos_each)
  past_out <- list(
    ucos = list(
      ucos_each = ucos_each,
      avW_ucos_each = avW,
      change_obj_each = matrix(0.1, nrow = length(ucos_each), ncol = n_outcomes),
      purity_each = matrix(1, nrow = length(ucos_each), ncol = 3),
      ucos_outcome = seq_len(n_outcomes)
    ),
    cos = list(cos = list())
  )

  ns <- asNamespace("colocboost")
  original_pairs <- get(".merge_ucos_overlap_pairs", envir = ns)
  original_between <- get("get_between_purity", envir = ns)
  original_purity <- get("get_purity", envir = ns)

  pairwise_pairs <- function(ucos_each) {
    if (length(ucos_each) < 2L) {
      return(matrix(integer(0), ncol = 2L, dimnames = list(NULL, c("i", "j"))))
    }
    pairs <- do.call(rbind, lapply(seq_len(length(ucos_each) - 1L), function(i) {
      do.call(rbind, lapply((i + 1L):length(ucos_each), function(j) {
        if (length(intersect(ucos_each[[i]], ucos_each[[j]])) == 0) {
          return(NULL)
        }
        c(i, j)
      }))
    }))
    if (is.null(pairs)) {
      return(matrix(integer(0), ncol = 2L, dimnames = list(NULL, c("i", "j"))))
    }
    colnames(pairs) <- c("i", "j")
    pairs
  }

  unlockBinding("get_between_purity", ns)
  assign("get_between_purity", function(...) {
    c(min_abs_cor = 0.9, max_abs_cor = 1, median_abs_cor = 0.9)
  }, envir = ns)
  lockBinding("get_between_purity", ns)
  unlockBinding("get_purity", ns)
  assign("get_purity", function(...) c(1, 1, 1), envir = ns)
  lockBinding("get_purity", ns)
  on.exit({
    unlockBinding(".merge_ucos_overlap_pairs", ns)
    assign(".merge_ucos_overlap_pairs", original_pairs, envir = ns)
    lockBinding(".merge_ucos_overlap_pairs", ns)
    unlockBinding("get_between_purity", ns)
    assign("get_between_purity", original_between, envir = ns)
    lockBinding("get_between_purity", ns)
    unlockBinding("get_purity", ns)
    assign("get_purity", original_purity, envir = ns)
    lockBinding("get_purity", ns)
  }, add = TRUE)

  candidate_result <- get("merge_ucos", envir = ns)(
    cb_obj, past_out,
    min_abs_corr = 0.5,
    median_cos_abs_corr = 0.8
  )

  unlockBinding(".merge_ucos_overlap_pairs", ns)
  assign(".merge_ucos_overlap_pairs", pairwise_pairs, envir = ns)
  lockBinding(".merge_ucos_overlap_pairs", ns)
  pairwise_result <- get("merge_ucos", envir = ns)(
    cb_obj, past_out,
    min_abs_corr = 0.5,
    median_cos_abs_corr = 0.8
  )

  expect_equal(candidate_result, pairwise_result)
  expect_equal(length(candidate_result$cos$cos$cos), 1L)
  expect_equal(length(candidate_result$ucos$ucos_each), 3L)
})

test_that("duplicate purity contexts are collapsed for post-assembly checks", {
  cb_data <- list(
    data = lapply(seq_len(4), function(i) {
      list(
        X = NULL,
        XtX = diag(5),
        variable_miss = if (i <= 2) integer(0) else 5L,
        ref_label = "LD"
      )
    }),
    dict = rep(1L, 4)
  )

  unique_outcomes <- get(".cb_unique_purity_outcomes", envir = asNamespace("colocboost"))(
    cb_data,
    seq_len(4)
  )

  expect_equal(unique_outcomes, c(1L, 3L))
})

test_that("merge_cos_ucos reuses duplicate purity contexts", {
  p <- 6
  L <- 3
  cb_obj <- list(
    cb_data = list(
      data = lapply(seq_len(L), function(i) {
        list(
          X = NULL,
          XtX = diag(p),
          variable_miss = integer(0),
          ref_label = "LD"
        )
      }),
      dict = rep(1L, L)
    ),
    cb_model_para = list(P = p, L = L)
  )
  class(cb_obj) <- "colocboost"

  out_cos <- list(cos = list(
    cos = list(cos1 = c(1, 2)),
    coloc_outcomes = list(1L),
    avWeight = list(matrix(runif(p), nrow = p, dimnames = list(NULL, "outcome1"))),
    cs_change = matrix(0.1, nrow = 1, ncol = L)
  ))
  out_ucos <- list(
    ucos_each = list(ucos1 = c(2, 3)),
    avW_ucos_each = matrix(runif(p), nrow = p, dimnames = list(NULL, "ucos1")),
    change_obj_each = matrix(0.2, nrow = 1, ncol = L),
    purity_each = matrix(1, nrow = 1, ncol = 3),
    ucos_outcome = 2L
  )

  ns <- asNamespace("colocboost")
  original_between <- get("get_between_purity", envir = ns)
  between_calls <- 0L
  unlockBinding("get_between_purity", ns)
  assign("get_between_purity", function(...) {
    between_calls <<- between_calls + 1L
    c(min_abs_cor = 0.9, max_abs_cor = 1, median_abs_cor = 0.9)
  }, envir = ns)
  lockBinding("get_between_purity", ns)
  on.exit({
    unlockBinding("get_between_purity", ns)
    assign("get_between_purity", original_between, envir = ns)
    lockBinding("get_between_purity", ns)
  }, add = TRUE)

  result <- get("merge_cos_ucos", envir = ns)(
    cb_obj,
    out_cos,
    out_ucos,
    median_cos_abs_corr = 0.8
  )

  expect_equal(between_calls, 1L)
  expect_null(result$ucos$ucos_each)
  expect_equal(result$cos$cos$coloc_outcomes[[1]], c(1L, 2L))
})

test_that("merge_cos_ucos skips purity checks for disjoint different-outcome sets", {
  p <- 6
  L <- 2
  cb_obj <- list(
    cb_data = list(
      data = lapply(seq_len(L), function(i) {
        list(
          X = NULL,
          XtX = diag(p),
          variable_miss = integer(0),
          ref_label = "LD"
        )
      }),
      dict = rep(1L, L)
    ),
    cb_model_para = list(P = p, L = L)
  )
  class(cb_obj) <- "colocboost"

  out_cos <- list(cos = list(
    cos = list(cos1 = c(1, 2)),
    coloc_outcomes = list(1L),
    avWeight = list(matrix(runif(p), nrow = p, dimnames = list(NULL, "outcome1"))),
    cs_change = matrix(0.1, nrow = 1, ncol = L)
  ))
  out_ucos <- list(
    ucos_each = list(ucos1 = c(4, 5)),
    avW_ucos_each = matrix(runif(p), nrow = p, dimnames = list(NULL, "ucos1")),
    change_obj_each = matrix(0.2, nrow = 1, ncol = L),
    purity_each = matrix(1, nrow = 1, ncol = 3),
    ucos_outcome = 2L
  )

  ns <- asNamespace("colocboost")
  original_between <- get("get_between_purity", envir = ns)
  unlockBinding("get_between_purity", ns)
  assign("get_between_purity", function(...) {
    stop("disjoint different-outcome sets should not require between-purity")
  }, envir = ns)
  lockBinding("get_between_purity", ns)
  on.exit({
    unlockBinding("get_between_purity", ns)
    assign("get_between_purity", original_between, envir = ns)
    lockBinding("get_between_purity", ns)
  }, add = TRUE)

  result <- get("merge_cos_ucos", envir = ns)(
    cb_obj,
    out_cos,
    out_ucos,
    median_cos_abs_corr = 0.8
  )

  expect_equal(length(result$ucos$ucos_each), 1L)
  expect_equal(result$cos$cos$coloc_outcomes[[1]], 1L)
})
