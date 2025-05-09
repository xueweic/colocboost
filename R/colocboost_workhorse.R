#' Workhorse of colocboost
#'
#' @details
#' The following functions are included in the gradient boosting iterations:
#'
#' Step 1: `colocboost_check_update_jk` selection scheme to choose outcomes need to be updated in each iteration.
#'
#' Step 2: `colocboost_update` core function of the proximity smoothed gradient boosting update.
#'
#' Step 3: check the stop criteria for each outcome.
#'
#' There is a version for LD free version with one causal assumption, implemented in `colocboost_one_causal`.
#'
#' @noRd
colocboost_workhorse <- function(cb_data,
                                 M = 500,
                                 tau = 0.01,
                                 prioritize_jkstar = TRUE,
                                 learning_rate_init = 0.01,
                                 learning_rate_decay = 1,
                                 func_simplex = "LD_z2z",
                                 lambda = 0.5,
                                 lambda_focal_outcome = 1,
                                 jk_equiv_corr = 0.8,
                                 jk_equiv_loglik = 1,
                                 stop_thresh = 1e-06,
                                 func_multi_test = "lfdr",
                                 stop_null = 0.05,
                                 multi_test_max = 1,
                                 multi_test_thresh = 1,
                                 ash_prior = "normal",
                                 p.adjust.methods = "fdr",
                                 func_compare = "min_max",
                                 coloc_thresh = 0.1,
                                 LD_free = FALSE,
                                 dynamic_learning_rate = TRUE,
                                 focal_outcome_idx = NULL,
                                 outcome_names = NULL) {
  if (!inherits(cb_data, "colocboost")) {
    stop("Input must from colocboost function!")
  }

  cb_model <- colocboost_init_model(cb_data,
    learning_rate_init = learning_rate_init,
    stop_thresh = stop_thresh,
    func_multi_test = func_multi_test,
    stop_null = stop_null,
    multi_test_max = multi_test_max,
    ash_prior = ash_prior,
    p.adjust.methods = p.adjust.methods,
    focal_outcome_idx = focal_outcome_idx
  )

  cb_model_para <- colocboost_init_para(cb_data, cb_model,
    tau = tau,
    func_simplex = func_simplex,
    lambda = lambda,
    lambda_focal_outcome = lambda_focal_outcome,
    learning_rate_decay = learning_rate_decay,
    multi_test_thresh = multi_test_thresh,
    multi_test_max = multi_test_max,
    func_multi_test = func_multi_test,
    LD_free = LD_free,
    outcome_names = outcome_names,
    focal_outcome_idx = focal_outcome_idx,
    dynamic_learning_rate = dynamic_learning_rate,
    prioritize_jkstar = prioritize_jkstar,
    jk_equiv_corr = jk_equiv_corr,
    jk_equiv_loglik = jk_equiv_loglik,
    func_compare = func_compare,
    coloc_thresh = coloc_thresh
  )
  
  if (M == 1){
    cb_model_para$M <- 1
  } else {
    M_single_outcome <- M
    M <- cb_model_para$L * M
    cb_model_para$M <- M
  }
  # - if multi_test value > multi_test cutoff for some outcomes, we will not update them.
  if (!is.null(cb_model_para$true_stop)) {
    if (sum(cb_model_para$update_y == 1) == 0) {
      # - if all outcomes do not have signals, STOP
      message(paste0(
        "Using multiple testing correction method: ", func_multi_test,
        ". Stop ColocBoost since no outcomes have association signals."
      ))
    } else {
      message(paste0(
        "Using multiple testing correction method: ", func_multi_test,
        ". Outcome ", paste(cb_model_para$true_stop, collapse = ", "),
        " for all variants are greater than ",
        multi_test_thresh, ". Will not update it!"
      ))
    }
    if (!is.null(focal_outcome_idx) & (M != 1)) {
      if (sum(cb_model_para$true_stop == focal_outcome_idx) != 0) {
        warning(paste("Stop ColocBoost since the focal outcome", focal_outcome_idx, "do not have association signals."))
        cb_model_para$update_y <- 0
      }
    }
  }

  if (sum(cb_model_para$update_y == 1) == 0) {
    cb_obj <- list("cb_data" = cb_data, "cb_model" = cb_model, "cb_model_para" = cb_model_para)
    class(cb_obj) <- "colocboost"
    return(cb_obj)
  }
  if (M == 1) {
    # single effect with or without LD matrix
    message("Running ColocBoost with assumption of one causal per outcome per region!")
    cb_obj <- colocboost_one_causal(cb_model, cb_model_para, cb_data)
    cb_obj$cb_model_para$coveraged <- "one_causal"
  } else {
    # - add more iterations for more outcomes
    for (m in 1:M) {
      if (sum(cb_model_para$update_y == 1) == 0) {
        break
      } else {
        # step 1: check which outcomes need to be updated at which variant
        cb_model_para <- colocboost_check_update_jk(cb_model, cb_model_para, cb_data)

        # step 2: gradient boosting for the updated outcomes
        cb_model <- colocboost_update(cb_model, cb_model_para, cb_data)

        # step 3: check stop for the updated ones
        # # - update cb_model and cb_model_parameter
        cb_model <- cb_model_update(cb_data, cb_model, cb_model_para)
        cb_model_para <- cb_model_para_update(cb_model, cb_model_para)

        # - stop by rtr < 0 - must stop
        pos.update <- which(cb_model_para$update_temp$update_status != 0)
        sum_cor <- sapply(cb_model[pos.update], function(cc) sum(cc$correlation))
        pos_rtr_stop <- which(sum_cor == 0)
        if (length(pos_rtr_stop) != 0) {
          cb_model_para$update_y[pos.update[pos_rtr_stop]] <- 0
          if (!is.null(focal_outcome_idx)){
            message_focal_text <- if (focal_outcome_idx %in% pos.update[pos_rtr_stop]) "including focal outcome" else NULL
          } else {message_focal_text <- NULL}
          message(paste(
            "Gradient boosting for outcome", paste(pos.update[pos_rtr_stop], collapse = ", "), message_focal_text,
            "stop since rtr < 0 or max(correlation) > 1 after", m, "iterations!",
            "Results for this locus are not stable, please check if mismatch between sumstat and LD!",
            "See details in tutorial website https://statfungen.github.io/colocboost/articles/."
          ))
        }

        # - stop by others - need to check iterations
        stop <- stop_no_coverage <- rep(NA, cb_model_para$L)
        for (i in pos.update) {
          if (i %in% pos.update[pos_rtr_stop]) {
            stop[i] <- stop_no_coverage[i] <- TRUE
            cb_model_para$coveraged_outcome[i] <- FALSE
          } else {
            stop_no_coverage[i] <- FALSE
            M_i <- length(cb_model[[i]]$profile_loglike_each)
            stop1 <- abs(cb_model[[i]]$profile_loglike_each[M_i] - cb_model[[i]]$profile_loglike_each[M_i - 1]) /
              cb_model[[i]]$profile_loglike_each[M_i - 1] < cb_model[[i]]$stop_thresh
            # stop2 = tail(abs(diff(cb_model[[i]]$obj_path)), n=1) < cb_model[[i]]$stop_thresh
            multiple_testing_correction <- get_multiple_testing_correction(
              z = cb_model[[i]]$z, miss_idx = cb_data$data[[i]]$variable_miss,
              func_multi_test = func_multi_test,
              ash_prior = ash_prior,
              p.adjust.methods = p.adjust.methods
            )
            cb_model[[i]]$multi_correction <- multiple_testing_correction
            stop2 <- all(multiple_testing_correction >= cb_model[[i]]$stop_null)
            stop3 <- (M_i > M_single_outcome)
            # -- to ensure if some outcomes do not update previously
            if (length(cb_model[[i]]$profile_loglike_each) >= 2) {
              stop[i] <- (stop1 | stop2 | stop3)
              if (stop3){
                stop_no_coverage[i] <- TRUE
                cb_model_para$coveraged_outcome[i] <- FALSE
                warning(paste("ColocBoost gradient boosting for outcome", i, "did not coverage in",
                              M_single_outcome, "iterations! Please check consistency between summary statistics",
                              "and LD matrix. See details in tutorial website https://statfungen.github.io/colocboost/articles/."))
                
              }
            } else {
              stop[i] <- FALSE
            }
          }
        }

        if (all(length(stop) == 1 & stop)) {
          cb_model_para$update_y <- 0
          if (cb_model_para$L == 1) {
            if (!stop_no_coverage) message(paste("Gradient boosting for outcome 1 converged after", m, "iterations!"))
          }
        } else {
          if (all(!stop[pos.update])) {
            # -- all Y do not stop, need more iterations.
            cb_model_para$update_y <- cb_model_para$update_y
          } else {
            pos_stop <- which(stop) # which outcome reach the stop criterion
            ttmp <- boost_check_stop(cb_model, cb_model_para, pos_stop, stop_no_coverage)
            cb_model_para <- ttmp$cb_model_para
            cb_model <- ttmp$cb_model
            # - if there is some outcomes need stop
            if (!is.null(cb_model_para$true_stop)) {
              ####### ---------------------------------------------
              # calculate objective function of Y for the last iteration.
              cb_model <- boost_obj_last(cb_data, cb_model, cb_model_para)
              tt_stop <- setdiff(cb_model_para$true_stop, cb_model_para$no_coverage_stop)
              if (!is.null(focal_outcome_idx)) {
                if (focal_outcome_idx %in% tt_stop) {
                  message(paste(
                    "Gradient boosting for focal outcome", focal_outcome_idx,
                    "converged after", m, "iterations!"
                  ))
                  if (length(setdiff(tt_stop, focal_outcome_idx)) != 0) {
                    message(paste(
                      "Gradient boosting for outcome", paste(setdiff(tt_stop, focal_outcome_idx), collapse = ", "),
                      "converged after", m, "iterations!"
                    ))
                  }
                } else {
                  if (length(tt_stop) != 0)
                    message(paste("Gradient boosting for outcome", paste(tt_stop, collapse = ", "), "converged after", m, "iterations!"))
                }
              } else {
                if (length(tt_stop) != 0)
                  message(paste("Gradient boosting for outcome", paste(tt_stop, collapse = ", "), "converged after", m, "iterations!"))
              }
            }
          }
        }
      }
      if (m %% 1000 == 0) {
        message(paste("Gradient boosting at", m, "iterations, still updating."))
      }
    }

    cb_model_para$num_updates <- m
    for (i in 1:length(cb_model)) {
      cb_model[[i]]$obj_path <- as.numeric(unlist(cb_model[[i]]$obj_path[-1]))
    }
    for (i in 1:length(cb_model)) {
      cb_model[[i]]$obj_single <- as.numeric(unlist(cb_model[[i]]$obj_single[-1]))
    }

    if (m == M) {
      ####### ---------------------------------------------
      # calculate objective function of Y for the last iteration.
      cb_model_para$true_stop <- which(cb_model_para$update_y == 1)
      ####### ---------------------------------------------
      # calculate objective function of Y for the last iteration.
      cb_model <- boost_obj_last(cb_data, cb_model, cb_model_para)
      warning(paste("ColocBoost updates did not converge in", M, "iterations; checkpoint at last iteration"))
      cb_model_para$coveraged <- FALSE
    }


    # -- remove redundant parameters
    rm_elements <- c("update_temp", "update_y")
    if (!is.null(cb_model_para$need_more)) {
      rm_elements <- c(rm_elements, "need_more")
    }
    if (!is.null(cb_model_para$true_stop)) {
      rm_elements <- c(rm_elements, "true_stop")
    }
    cb_model_para[rm_elements] <- NULL

    cb_obj <- list("cb_data" = cb_data, "cb_model" = cb_model, "cb_model_para" = cb_model_para)
    class(cb_obj) <- "colocboost"
  }
  return(cb_obj)
}

# - update cb_model at the begining
cb_model_update <- function(cb_data, cb_model, cb_model_para) {
  pos <- which(cb_model_para$update_temp$update_status != 0)
  for (i in pos) {
    X_dict <- cb_data$dict[i]
    data_each <- cb_data$data[[i]]
    model_each <- cb_model[[i]]
    tmp <- get_correlation(
      X = cb_data$data[[X_dict]]$X, res = model_each$res, XtY = data_each$XtY,
      N = data_each$N, YtY = data_each$YtY,
      XtX = cb_data$data[[X_dict]]$XtX,
      beta_k = model_each$beta,
      miss_idx = data_each$variable_miss
    )
    cb_model[[i]]$correlation <- tmp
    cb_model[[i]]$z <- get_z(tmp, n = data_each$N, model_each$res)
  }
  return(cb_model)
}


cb_model_para_update <- function(cb_model, cb_model_para) {
  # - additional calculating profile_loglikelihood
  tmp <- sum(sapply(1:length(cb_model), function(i) tail(cb_model[[i]]$profile_loglike_each, n = 1)))
  cb_model_para$profile_loglike <- c(cb_model_para$profile_loglike, tmp)
  num_updates_outcomes <- sapply(1:length(cb_model), function(i) length(cb_model[[i]]$profile_loglike_each))
  names(num_updates_outcomes) <- names(cb_model)
  cb_model_para$num_updates_outcome <- num_updates_outcomes
  return(cb_model_para)
}
