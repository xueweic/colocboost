#!/usr/bin/env Rscript
# Benchmark: Phase 1 optimization (XtX*beta cache + precomputed constants)
#
# This script compares the optimized code against a simulated "no-cache" baseline
# by manually running the dominant O(P^2) operations the number of times they
# would occur with and without caching.
#
# The optimization eliminates 2 of 3 redundant XtX %*% beta computations per
# iteration per outcome.

library(MASS)

cat("=== ColocBoost Phase 1 Optimization Benchmark ===\n\n")

# ---- Generate test data at different scales ----
run_benchmark <- function(p, L, M, n_ref = 500, seed = 42) {
  set.seed(seed)

  cat(sprintf("P = %d variants, L = %d outcomes, M = %d iterations, N_ref = %d\n", p, L, M, n_ref))

  # Generate LD matrix (P x P)
  sigma <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigma[i, j] <- 0.9^abs(i - j)
    }
  }
  # Ensure positive definite
  LD <- sigma

  # Generate random beta vector
  beta <- rnorm(p) * 0.01

  cat(sprintf("  LD matrix size: %.1f MB\n", object.size(LD) / 1024^2))

  # ---- Benchmark: XtX %*% beta ----
  # Before optimization: 3 * M * L calls to XtX %*% beta
  # After optimization: 1 * M * L calls to XtX %*% beta

  n_calls_before <- 3 * M * L
  n_calls_after <- 1 * M * L

  # Time a single XtX %*% beta
  t_single <- system.time({
    for (rep in 1:100) {
      result <- LD %*% beta
    }
  })[["elapsed"]] / 100

  cat(sprintf("  Single XtX %%*%% beta: %.4f seconds\n", t_single))
  cat(sprintf("  Before optimization: %d calls = %.2f seconds\n",
              n_calls_before, n_calls_before * t_single))
  cat(sprintf("  After optimization:  %d calls = %.2f seconds\n",
              n_calls_after, n_calls_after * t_single))
  cat(sprintf("  Speedup on dominant cost: %.1fx\n",
              n_calls_before * t_single / (n_calls_after * t_single)))
  cat(sprintf("  Time saved: %.2f seconds\n\n",
              (n_calls_before - n_calls_after) * t_single))

  invisible(list(
    p = p, L = L, M = M,
    t_single = t_single,
    t_before = n_calls_before * t_single,
    t_after = n_calls_after * t_single
  ))
}

# ---- Run end-to-end colocboost benchmark ----
run_colocboost_benchmark <- function(p, L, M, n = 200, seed = 42) {
  set.seed(seed)

  cat(sprintf("\n--- End-to-end colocboost: P=%d, L=%d, M=%d ---\n", p, L, M))

  sigma <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigma[i, j] <- 0.9^abs(i - j)
    }
  }
  X <- mvrnorm(n, rep(0, p), sigma)
  colnames(X) <- paste0("SNP", 1:p)
  Y <- matrix(rnorm(n * L), n, L)
  # Add signal to SNP5 for all traits
  for (l in 1:L) {
    Y[, l] <- Y[, l] + X[, 5] * 0.5
  }
  LD <- cor(X)

  # Generate summary statistics
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
      beta = beta_hat, sebeta = se_hat, z = z,
      n = n, variant = colnames(X)
    )
  }

  # Time full colocboost run
  t_full <- system.time({
    suppressWarnings(suppressMessages({
      result <- colocboost::colocboost(
        sumstat = sumstat_list,
        LD = LD,
        M = M,
        output_level = 1
      )
    }))
  })[["elapsed"]]

  cat(sprintf("  Total wall time: %.2f seconds\n", t_full))
  invisible(t_full)
}

# ---- Scenarios ----

cat("--- Micro-benchmark: XtX %*% beta operation ---\n\n")

results <- list()
results[[1]] <- run_benchmark(p = 1000, L = 2, M = 100)
results[[2]] <- run_benchmark(p = 2000, L = 5, M = 200)
results[[3]] <- run_benchmark(p = 5000, L = 3, M = 300)
results[[4]] <- run_benchmark(p = 5000, L = 10, M = 500)

cat("\n--- Summary Table ---\n")
cat(sprintf("%-8s %-4s %-5s %-12s %-12s %-8s\n",
            "P", "L", "M", "Before(s)", "After(s)", "Speedup"))
for (r in results) {
  cat(sprintf("%-8d %-4d %-5d %-12.2f %-12.2f %-8.1fx\n",
              r$p, r$L, r$M, r$t_before, r$t_after, r$t_before / r$t_after))
}

cat("\n--- End-to-end colocboost timings ---\n")
run_colocboost_benchmark(p = 100, L = 2, M = 50)
run_colocboost_benchmark(p = 100, L = 5, M = 100)
