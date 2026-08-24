library(testthat)

test_that("colocalization candidate grouping supports more than 10000-byte signatures", {
  n_outcomes <- 4589L
  update <- matrix(-1L, nrow = n_outcomes, ncol = 3L)
  update[seq_len(100L), 1:2] <- 1L
  update[101:200, 3] <- 1L

  old_signature <- paste0(update[, 1], collapse = ",")
  expect_gt(nchar(old_signature, type = "bytes"), 10000L)

  result <- .group_coloc_candidates(update, seq_len(ncol(update)))

  expect_equal(unname(result), list(c(1L, 2L), 3L))
  expect_true(all(nchar(names(result), type = "bytes") < 10000L))
})

test_that("CoS details support more than 10000-byte outcome-set names", {
  n_outcomes <- 4589L
  cos_name <- paste0("cos1:", paste0("y", seq_len(n_outcomes), collapse = "_"))
  expect_gt(nchar(cos_name, type = "bytes"), 10000L)

  local_mocked_bindings(
    get_cos_evidence = function(...) {
      list(
        normalization_evidence = list(data.frame(npc_outcome = 1)),
        npc = 1
      )
    },
    get_integrated_weight = function(...) c(1, 0),
    get_in_cos = function(...) list(1L),
    get_purity = function(...) c(1, 1, 1),
    .package = "colocboost"
  )

  data_entry <- list(
    X = matrix(0, nrow = 2, ncol = 2),
    XtX = NULL,
    N = 2L,
    variable_miss = integer(0),
    ref_label = "individual"
  )
  cb_obj <- list(
    cb_model_para = list(
      outcome_names = paste0("trait", seq_len(n_outcomes)),
      variables = c("v1", "v2"),
      P = 2L,
      coverage = 0.95,
      min_abs_corr = 0.5,
      median_abs_corr = NULL,
      n_purity = 100L,
      weight_fudge_factor = 1.5,
      use_entropy = FALSE,
      residual_correlation = NULL
    ),
    cb_data = list(
      dict = rep(1L, n_outcomes),
      data = rep(list(data_entry), n_outcomes)
    )
  )
  coloc_out <- list(
    cos = list(cos1 = 1L),
    coloc_outcomes = list(seq_len(n_outcomes)),
    avWeight = list(matrix(1, nrow = 2L, ncol = n_outcomes)),
    cs_change = matrix(1, nrow = 1L, ncol = n_outcomes)
  )
  data_info <- list(variables = c("v1", "v2"))

  expect_no_error(result <- get_cos_details(cb_obj, coloc_out, data_info))
  expect_equal(names(result$cos_results$cos$cos_index), cos_name)
})

test_that("robust CoS filtering supports more than 10000-byte outcome-set names", {
  n_outcomes <- 4589L
  outcome_idx <- seq_len(n_outcomes)
  outcome_names <- paste0("trait", outcome_idx)
  weights <- matrix(rep(c(0.9, 0.1), n_outcomes), nrow = 2L)
  colnames(weights) <- paste0("outcome", outcome_idx)
  outcome_npc <- data.frame(
    relative_logLR = rep(1, n_outcomes),
    npc_outcome = rep(1, n_outcomes),
    outcomes_index = outcome_idx,
    row.names = outcome_names
  )
  purity <- matrix(
    1,
    2L,
    2L,
    dimnames = list(c("cos1", "cos2"), c("cos1", "cos2"))
  )

  cb_output <- list(
    data_info = list(
      variables = c("v1", "v2"),
      z = rep(list(c(10, 0)), n_outcomes),
      coef = rep(list(c(1, 0)), n_outcomes),
      n_outcomes = n_outcomes,
      n_variables = 2L,
      outcome_info = data.frame(
        outcome_names = outcome_names,
        is_focal = rep(FALSE, n_outcomes)
      )
    ),
    vcp = c(v1 = 0.9, v2 = 0.1),
    cos_details = list(
      cos_top_variables = data.frame(
        top_index = c(1L, 2L),
        top_variables = c("v1", "v2"),
        row.names = c("cos1", "cos2")
      ),
      cos = list(
        cos_index = list(cos1 = 1L, cos2 = 2L),
        cos_variables = list(cos1 = "v1", cos2 = "v2")
      ),
      cos_outcomes = list(
        outcome_index = list(cos1 = outcome_idx, cos2 = 1L),
        outcome_name = list(cos1 = outcome_names, cos2 = outcome_names[1])
      ),
      cos_outcomes_npc = list(
        cos1 = outcome_npc,
        cos2 = data.frame(
          relative_logLR = 0,
          npc_outcome = 0,
          outcomes_index = 1L,
          row.names = outcome_names[1]
        )
      ),
      cos_vcp = list(cos1 = c(0.9, 0.1), cos2 = c(0.1, 0.9)),
      cos_weights = list(
        cos1 = weights,
        cos2 = matrix(c(0.1, 0.9), ncol = 1L, dimnames = list(NULL, "outcome1"))
      ),
      cos_npc = c(cos1 = 1, cos2 = 1),
      cos_min_npc_outcome = c(cos1 = 1, cos2 = 1),
      cos_purity = list(
        min_abs_cor = purity,
        median_abs_cor = purity,
        max_abs_cor = purity
      )
    )
  )
  class(cb_output) <- "colocboost"

  expect_no_error(
    result <- get_robust_colocalization(
      cb_output,
      cos_npc_cutoff = 0,
      npc_outcome_cutoff = 0,
      pvalue_cutoff = 1
    )
  )
  expect_equal(length(result$cos_details$cos_outcomes$outcome_index[[1]]), n_outcomes)
  expect_gt(nchar(names(result$cos_details$cos$cos_index)[1], type = "bytes"), 10000L)
  expect_no_error(plot_input <- get_input_plot(result))
  expect_equal(length(plot_input$cos_vcp), n_outcomes)
})
