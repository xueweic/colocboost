#' @title Set of functions for post inferences from ColocBoost
#'
#' @description
#' The `colocboost_get_methods` functions access basic properties inferences from a fitted ColocBoost model. This documentation serves as a summary for all related post-inference functions.
#'
#'
#' @details
#' The following functions are included in this set:
#' `get_summary_table` get the colocalization summary table with or without the specific outcomes.
#'
#' These functions are not exported individually and are accessed via `colocboost_get_methods`.
#'
#' @rdname colocboost_get_methods
#' @keywords cb_post_inference
#' @export
colocboost_get_methods <- function() {
    message("This function post inferences of colocboost output. See details for more information.")
}

#' @noRd
#' @keywords cb_post_inference
get_data_info <- function(cb_obj){

    ## - analysis data information
    n_outcome <- cb_obj$cb_model_para$L
    n_variables <- cb_obj$cb_model_para$P
    analysis_outcome <- cb_obj$cb_model_para$outcome_names
    variables <- cb_obj$cb_data$variable.names
    target_outcome <- NULL
    is_target <- rep(FALSE, n_outcome)
    if (!is.null(cb_obj$cb_model_para$target_idx)){
      target_outcome <- analysis_outcome[cb_obj$cb_model_para$target_idx]
      is_target[cb_obj$cb_model_para$target_idx] <- TRUE
    }
    is_sumstat <- grepl("sumstat_outcome", names(cb_obj$cb_data$data))
    outcome_info <- data.frame(
      "outcome_names" = analysis_outcome,  "sample_size" = cb_obj$cb_model_para$N,
      "is_sumstats" = is_sumstat, "is_target" = is_target )
    rownames(outcome_info) <- paste0("y", 1:n_outcome)
    
    ## - marginal associations
    z_scores <- lapply(cb_obj$cb_model, function(cb){ as.numeric(cb$z_univariate) })
    betas <- lapply(cb_obj$cb_model, function(cb){ as.numeric(cb$beta) })
    names(z_scores) <- names(betas) <- analysis_outcome
    if (all(grepl("chr", variables))){
       # need to check with AI
       # chr:pos:a1:a2
       # 1:pos:a1
       # chr_pos_
       # chr.pos
       # y1: A, B, D, F, G
       # y2: B, C, D, E
       # A,B,C,D,E,F,G
      
        # - if variables_name has position informtaion, re-order all_info based on positions
        # position <- as.numeric(sapply(variables, function(tmp) strsplit(tmp, ":")[[1]][2]))
        position <- as.numeric(gsub(".*:(\\d+).*", "\\1", variables))
        ordered <- order(position, decreasing = FALSE)
        z_scores <- lapply(z_scores, function(z) z[ordered])
        variables <- variables[ordered]
    }

    ## - output data info
    data.info <- list("n_outcomes" = n_outcome,
                      "n_variables" = n_variables,
                      "outcome_info" = outcome_info,
                      "variables" = variables,
                      "coef" = betas,
                      "z" = z_scores)
    return(data.info)

}



#' @noRd
#' @keywords cb_post_inference
get_cos_details <- function(cb_obj, coloc_out, data_info = NULL, npc_cutoff = 0.7){

    if (is.null(data_info))
        data_info <- get_data_info(cb_obj)


    ### ----- Define the colocalization results
    coloc_sets <- coloc_out$cos
    if (length(coloc_sets)!=0){

      # - colocalization outcome configurations
      tmp <- get_cos_evidence(cb_obj, coloc_out, data_info, npc_cutoff = npc_cutoff)
      normalization_evidence <- tmp$normalization_evidence
      npc <- tmp$npc
      
      # - colocalized outcomes
        analysis_outcome <- cb_obj$cb_model_para$outcome_names
        coloc_outcome_index <- coloc_outcome <- list()
        colocset_names <- c()
        for (i in 1:length(coloc_out$cos)){
            coloc_outcome_index[[i]] <- coloc_out$coloc_outcomes[[i]]
            coloc_outcome[[i]] <- analysis_outcome[coloc_outcome_index[[i]]]
            colocset_names[i] <- paste0("cos", i, ":", paste0(paste0("y", coloc_outcome_index[[i]]), collapse = "_"))
            if (grepl("merged", names(coloc_sets)[i])) {
              colocset_names[i] <- paste0(colocset_names[i], ":merged")
            }
        }
        names(coloc_outcome) <- names(coloc_outcome_index) <- colocset_names
        names(npc) <- names(normalization_evidence) <- colocset_names
        coloc_outcomes <- list("outcome_index" = coloc_outcome_index, "outcome_name" = coloc_outcome)

        # - colocalized sets for variables
        coloc_csets_variableidx <- coloc_out$cos
        coloc_csets_variablenames <- lapply(coloc_csets_variableidx, function(coloc_tmp){ cb_obj$cb_model_para$variables[coloc_tmp] })
        coloc_csets_variableidx <- lapply(coloc_csets_variablenames, function(variable) match(variable, data_info$variables))
        names(coloc_csets_variableidx) <- names(coloc_csets_variablenames) <- colocset_names
        coloc_csets_original <- list("cos_index" = coloc_csets_variableidx, "cos_variables" = coloc_csets_variablenames)

        # - colocalized set cs_change
        cs_change <- coloc_out$cs_change
        rownames(cs_change) <- colocset_names
        colnames(cs_change) <- analysis_outcome

        # - VCP
        cos_weights <- lapply(coloc_out$avWeight, function(w){
          pos <- match(data_info$variables, cb_obj$cb_model_para$variables)
          return(w[pos, , drop = FALSE])
        })
        int_weight <- lapply(cos_weights, get_integrated_weight, alpha = cb_obj$cb_model_para$alpha)
        names(int_weight) <- names(cos_weights) <- colocset_names
        vcp <- as.vector(1 - apply(1 - do.call(cbind, int_weight),1,prod))
        names(vcp) <- data_info$variables
        
        # - resummary results
        cos_re_idx <- lapply(int_weight, function(w){ unlist(get_in_csets(w, coverage = cb_obj$cb_model_para$coverage))})
        cos_re_var <- lapply(cos_re_idx, function(idx){ data_info$variables[idx] })
        coloc_csets <- list("cos_index" = cos_re_idx, "cos_variables" = cos_re_var)
        
        # - hits variables in each csets
        coloc_hits <- coloc_hits_variablenames <- coloc_hits_names <- c()
        for (i in 1:length(int_weight)){
            inw <- int_weight[[i]]
            pp <- which(inw == max(inw))
            coloc_hits <- c(coloc_hits, pp)
            coloc_hits_variablenames <- c(coloc_hits_variablenames, data_info$variables[pp])
            if (length(pp)==1){
                coloc_hits_names <- c(coloc_hits_names, names(int_weight)[i])
            } else {
                coloc_hits_names <- c(coloc_hits_names, paste0(names(int_weight)[i], ".", 1:length(pp)))
            }
        }
        coloc_hits <- data.frame("top_index" = coloc_hits, "top_variables" = coloc_hits_variablenames) 
        rownames(coloc_hits) <- coloc_hits_names

        # - purity
        ncos <- length(coloc_csets$cos_index)
        if (ncos >= 2){
            empty_matrix <- matrix(NA, ncos, ncos)
            colnames(empty_matrix) <- rownames(empty_matrix) <- colocset_names
            csets_purity <- lapply(1:3, function(ii){
              diag(empty_matrix) <- coloc_out$purity[,ii]
              return(empty_matrix)
            })
            for (i in 1:(ncos-1)){
                for (j in (i+1):ncos){
                    cset1 <- coloc_csets$cos_index[[i]]
                    cset2 <- coloc_csets$cos_index[[j]]
                    y.i <- coloc_outcomes$outcome_index[[i]]
                    y.j <- coloc_outcomes$outcome_index[[j]]
                    yy <- unique(y.i, y.j)
                    res <- list()
                    flag <- 1
                    for (ii in yy){
                        X_dict <- cb_obj$cb_data$dict[ii]
                        res[[flag]] <- get_between_purity(cset1, cset2, X = cb_obj$cb_data$data[[X_dict]]$X,
                                                          Xcorr = cb_obj$cb_data$data[[X_dict]]$XtX,
                                                          N = cb_obj$cb_data$data[[ii]]$N,
                                                          miss_idx = cb_obj$cb_data$data[[ii]]$variable_miss,
                                                          P = cb_obj$cb_model_para$P)
                        flag <- flag + 1
                    }
                    res <- Reduce(pmax, res)
                    csets_purity <- lapply(1:3, function(ii){
                        csets_purity[[ii]][i,j] <- csets_purity[[ii]][j,i] <- res[ii]
                        return(csets_purity[[ii]])
                    })
                }
            }
            names(csets_purity) <- c("min_abs_cor", "max_abs_cor", "median_abs_cor")

        } else {
            csets_purity <- coloc_out$purity
            rownames(csets_purity) <- colocset_names
        }

        # - save coloc_results
        coloc_results <- list("cos" = coloc_csets,
                              "cos_outcomes" = coloc_outcomes,
                              "cos_vcp" = int_weight,
                              "cos_npc" = npc,
                              "cos_purity" = csets_purity,
                              "cos_top_variables" = coloc_hits,
                              "cos_outcomes_delta" = cs_change,
                              "cos_normalization_evidence" = normalization_evidence,
                              "cos_weights" = cos_weights)
        
        
        # - missing variable and warning message
        missing_variables_idx <- Reduce(union, lapply(cb_obj$cb_data$data,function(cb) cb$variable_miss ))
        missing_variables <- cb_obj$cb_model_para$variables[missing_variables_idx]
        cos_missing_variables_idx <-  lapply(coloc_csets_original$cos_variables, function(variable){
          missing <- intersect(variable, missing_variables)
          if (length(missing)!=0){
            match(missing, data_info$variables)
          } else { NULL }
        })
        cos_missing_variables <- lapply(cos_missing_variables_idx, function(variable){if(!is.null(variable)){data_info$variables[variable]} else {NULL}})
        warning_needed <- any(!sapply(cos_missing_variables, is.null))
        if (warning_needed){
           is_missing <- which(!sapply(cos_missing_variables, is.null))
           cos_missing_variables_idx <- cos_missing_variables_idx[is_missing]
           cos_missing_variables <- cos_missing_variables[is_missing]
           cos_missing_vcp <- lapply(cos_missing_variables_idx, function(idx){ vcp[idx] })
           warning_message <- paste("CoS", paste(names(cos_missing_variables_idx), collapse = ","),
                                    "contains missing variables in at least one outcome.",
                                    "The missing variables will cause the ~0 VCP scores.")
           cos_warnings <- list("cos_missing_info" = list("index" = cos_missing_variables_idx,
                                                          "variables" = cos_missing_variables,
                                                          "vcp" = cos_missing_vcp),
                                "warning_message" = warning_message)
           coloc_results$cos_warnings <- cos_warnings
        }

    } else {
        coloc_results <- NULL
        vcp <- NULL
    }
    return(list("cos_results"=coloc_results, "vcp"=vcp))

}

#' @noRd
#' @keywords cb_post_inference
get_model_info <- function(cb_obj, outcome_names = NULL){
  
  if (is.null(outcome_names)){
    data_info <- get_data_info(cb_obj)
    outcome_names <- data_info$outcome_info$outcome_names
  }
  
  profile_loglik <- cb_obj$cb_model_para$profile_loglike
  n_updates <- cb_obj$cb_model_para$num_updates
  model_coveraged <- cb_obj$cb_model_para$coveraged
  jk_update <- cb_obj$cb_model_para$real_update_jk
  outcome_proximity_obj <- lapply(cb_obj$cb_model, function(cb) cb$obj_path)
  outcome_coupled_obj <- lapply(cb_obj$cb_model, function(cb) cb$obj_single)
  outcome_profile_loglik <- lapply(cb_obj$cb_model, function(cb) cb$profile_loglike_each)
  names(outcome_proximity_obj) <- names(outcome_coupled_obj) <-
    names(outcome_profile_loglik) <- outcome_names
  ll <- list("model_coveraged" = model_coveraged,
             "n_updates" = n_updates,
             "profile_loglik" = profile_loglik,
             "outcome_profile_loglik" = outcome_profile_loglik,
             "outcome_proximity_obj" = outcome_proximity_obj,
             "outcome_coupled_obj" = outcome_coupled_obj,
             "jk_update" = jk_update)
  return(ll)
  
}



#' @rdname colocboost_get_methods
#' @title Extract colocalization summary table
#'
#' @description
#' Get the colocalization summary table with or without the specific outcomes.
#'
#' @examples
#' get_summary_table(cb_output)
#' get_summary_table(cb_output, target_outcome = c("Y1", "Y2"))
#'
#' @noRd
#' @keywords cb_post_inference
get_cos_summary <- function(cb_output, outcome_names = NULL, gene_name = NULL, 
                            target_outcome = NULL){

    coloc_csets <- cb_output$cos_details$cos$cos_index
    if (length(coloc_csets) != 0){

        analysis_outcome <- cb_output$data_info$outcome_info$outcome_names
        if (!is.null(outcome_names)){ analysis_outcome <- outcome_names }
        coloc_outcome <- lapply(cb_output$cos_details$cos_outcomes$outcome_index, function(idx) analysis_outcome[idx])
        coloc_sets <- cb_output$cos_details$cos$cos_index
        if (!is.null(cb_output$cos_warnings)){
            cos_warnings
        }
        vcp <- as.numeric(cb_output$vcp)
    
        summary_table <- matrix(NA, nrow = length(coloc_sets), ncol = 11)
        colnames(summary_table) <- c("target_outcome", "colocalized_outcomes", "cos_id", "purity", 
                                     "top_variable", "top_variable_vcp", "cos_npc", "n_variables", 
                                     "colocalized_index", "colocalized_variables", "colocalized_variables_vcp")
        summary_table <- as.data.frame(summary_table)
        summary_table[,1] <- FALSE
        summary_table[,2] <- unlist(sapply(coloc_outcome, function(tmp) paste0(tmp, collapse = "; ")))
        summary_table[,3] <- names(coloc_sets)
        summary_table[,4] <- as.numeric(diag(as.matrix(cb_output$cos_details$cos_purity$min_abs_cor)))
        summary_table[,5] <- unlist(sapply(cb_output$cos_details$cos$cos_variables, function(tmp) tmp[1]))    
        summary_table[,6] <- sapply(coloc_sets, function(tmp) max(vcp[tmp]))   
        summary_table[,7] <- as.numeric(cb_output$cos_details$cos_npc)
        summary_table[,8] <- as.numeric(sapply(coloc_sets, length))
        summary_table[,9] <- unlist(sapply(coloc_sets, function(tmp) paste0(tmp, collapse = "; ")))
        summary_table[,10] <- unlist(sapply(cb_output$cos_details$cos$cos_variables, function(tmp) paste0(tmp, collapse = "; ")))
        summary_table[,11] <- unlist(sapply(coloc_sets, function(tmp) paste0(vcp[tmp], collapse = "; ")))
        if (!is.null(gene_name)){ summary_table$gene_name <- gene_name }
        
        if (!is.null(target_outcome)){
            tmp <- sapply(target_outcome, function(tmp) grep(paste0(tmp, "\\b"), analysis_outcome))
            if (all(sapply(tmp, length)!=0)){
                if.target <- sapply(coloc_outcome, function(cp){
                    tt <- sapply(target_outcome, function(tmp) grep(paste0(tmp, "\\b"), cp))
                    all(sapply(tt, length)!=0)
                })
                summary_table$target_outcome <- ifelse(if.target, target_outcome, FALSE)
                summary_table=summary_table[order(summary_table$target_outcome == "FALSE"), ]
                if (sum(if.target) == 0){ warnings("No colocalization with target outcomes.") }
            } else {
                warnings("Target outcome is not in the analysis outcomes, please check.")
            }
        }
    } else {
        summary_table <- NULL
    }
    return(summary_table)

}


#' @noRd
#' @keywords cb_post_inference
get_full_output <- function(cb_obj, past_out = NULL, variables = NULL, cb_output = NULL, weaker_ucos = FALSE){

    cb_model <- cb_obj$cb_model
    cb_model_para <- cb_obj$cb_model_para

    ## - obtain the order of variables based on the variables names if it has position information
    if (is.null(variables)){
        variables <- cb_obj$cb_model_para$variables
        if (all(grepl("chr", variables))){
            # - if variables_name has position informtaion, re-order all_info based on positions
            position <- as.numeric(gsub(".*:(\\d+).*", "\\1", variables))
            # position <- as.numeric(sapply(variables_name, function(tmp) strsplit(tmp, ":")[[1]][2]))
            ordered <- order(position, decreasing = FALSE)
        } else {
            ordered <- 1:length(variables)
        }
    } else {
        ordered <- match(variables, cb_obj$cb_model_para$variables)
    }

    ## - reorder all output
    # - cb_model
    tmp <- lapply(cb_model, function(cb){
        cb$beta <- cb$beta[ordered]
        cb$weights_path <- cb$weights_path[, ordered]
        cb$change_loglike <- cb$change_loglike[ordered]
        cb$correlation <- as.numeric(cb$correlation[ordered])
        cb$z <- as.numeric(cb$z[ordered])
        cb$ld_jk <- cb$ld_jk[,ordered]
        cb$z_univariate <- as.numeric(cb$z_univariate[ordered])
        cb$beta_hat <- as.numeric(cb$beta_hat[ordered])
        cb$multi_correction <- as.numeric(cb$multi_correction[ordered])
        cb$multi_correction_univariate <- as.numeric(cb$multi_correction_univariate[ordered])
        return(cb)
    })
    cb_model <- tmp

    # - sets
    if (!is.null(past_out)){
        out_ucos <- past_out$ucos
        # - single sets
        if (!is.null(out_ucos$ucos_each)){
        
            out_ucos$ucos_each <- lapply(out_ucos$ucos_each, function(cs){
                match(cb_model_para$variables[cs], variables)
            })
            out_ucos$avW_ucos_each <- out_ucos$avW_ucos_each[ordered,,drop=FALSE]
            
            # - re-orginize specific results
            analysis_outcome <- cb_obj$cb_model_para$outcome_names
            specific_outcome_index <- specific_outcome <- list()
            specific_cs_names <- c()
            for (i in 1:length(out_ucos$ucos_each)){
              cc <- out_ucos$avW_ucos_each[,i,drop=FALSE] 
              tmp_names <- colnames(cc)
              specific_outcome_index[[i]] <- as.numeric(gsub(".*Y([0-9]+).*", "\\1", tmp_names))
              specific_outcome[[i]] <- analysis_outcome[specific_outcome_index[[i]]]
              specific_cs_names[i] <- paste0("ucos", i, ":y", specific_outcome_index[[i]])
            }
            names(specific_outcome) <- names(specific_outcome_index) <- specific_cs_names
            specific_outcomes <- list("outcome_index" = specific_outcome_index, "outcome_name" = specific_outcome)
            
            # - specific sets for variables
            specific_cs_variableidx <- out_ucos$ucos_each
            specific_cs_variablenames <- lapply(specific_cs_variableidx, function(specific_tmp){ cb_obj$cb_model_para$variables[specific_tmp] })
            specific_cs_variableidx <- lapply(specific_cs_variablenames, function(variable) match(variable, variables))
            names(specific_cs_variableidx) <- names(specific_cs_variablenames) <- specific_cs_names
            specific_css <- list("ucos_index" = specific_cs_variableidx, "ucos_variables" = specific_cs_variablenames)
            
            # - specific set cs_change
            cs_change <- out_ucos$change_obj_each
            rownames(cs_change) <- specific_cs_names
            colnames(cs_change) <- analysis_outcome
            index_change <- as.data.frame(which(cs_change != 0, arr.ind = TRUE))
            change_outcomes <- analysis_outcome[index_change$col]
            change_values <- diag(as.matrix(cs_change[index_change$row, index_change$col]))
            cs_change <- data.frame("ucos_outcome" = change_outcomes, "ucos_delta" = change_values)
            
            # - filter weak ucos
            check_null_max <- sapply(cb_model, function(cb) cb$check_null_max)
            remove_weak <- sapply(1:nrow(cs_change), function(ic){
              outcome_tmp <- cs_change$ucos_outcome[ic]
              delta_tmp <- cs_change$ucos_delta[ic]
              pp <- which(cb_obj$cb_model_para$outcome_names == outcome_tmp)
              check_tmp <- check_null_max[pp]
              delta_tmp >= check_tmp
            })
            keep_ucos <- which(remove_weak)
            if (length(keep_ucos)==0){
              specific_results <- NULL
            } else {
              specific_outcomes$outcome_index <- specific_outcomes$outcome_index[keep_ucos]
              specific_outcomes$outcome_name <- specific_outcomes$outcome_name[keep_ucos]
              specific_css$ucos_index <- specific_css$ucos_index[keep_ucos]
              specific_css$ucos_variables <- specific_css$ucos_variables[keep_ucos]
              cs_change <- cs_change[keep_ucos, , drop = FALSE]
              out_ucos$avW_ucos_each <- out_ucos$avW_ucos_each[, keep_ucos, drop = FALSE]
              specific_cs_names <- specific_cs_names[keep_ucos]
              out_ucos$purity_each <- out_ucos$purity_each[keep_ucos, , drop = FALSE]
              
              # - ucos_weight
              specific_w <- lapply(1:ncol(out_ucos$avW_ucos_each), function(ii) out_ucos$avW_ucos_each[,ii,drop=FALSE])
              names(specific_w) <- specific_cs_names
              
              # - hits variables in each csets
              cs_hits <- sapply(1:length(specific_w), function(jj){ inw=specific_w[[jj]]; sample(which(inw == max(inw)),1) })
              cs_hits_variablenames <- sapply(cs_hits, function(ch) variables[ch] )
              specific_cs_hits <- data.frame("top_index" = cs_hits, "top_variables" = cs_hits_variablenames)   # save
              rownames(specific_cs_hits) <- specific_cs_names
              
              # - purity
              nucos <- length(specific_css$ucos_index)
              if (nucos >= 2){
                empty_matrix <- matrix(NA, nucos, nucos)
                colnames(empty_matrix) <- rownames(empty_matrix) <- specific_cs_names
                specific_cs_purity <- lapply(1:3, function(ii){
                  diag(empty_matrix) <- out_ucos$purity_each[,ii]
                  return(empty_matrix)
                })
                for (i in 1:(nucos-1)){
                  for (j in (i+1):nucos){
                    cset1 <- specific_css$ucos_index[[i]]
                    cset2 <- specific_css$ucos_index[[j]]
                    y.i <- specific_outcomes$outcome_index[[i]]
                    y.j <- specific_outcomes$outcome_index[[j]]
                    yy <- unique(y.i, y.j)
                    res <- list()
                    flag <- 1
                    for (ii in yy){
                      X_dict <- cb_obj$cb_data$dict[ii]
                      res[[flag]] <- get_between_purity(cset1, cset2, X = cb_obj$cb_data$data[[X_dict]]$X,
                                                        Xcorr = cb_obj$cb_data$data[[X_dict]]$XtX,
                                                        N = cb_obj$cb_data$data[[ii]]$N,
                                                        miss_idx = cb_obj$cb_data$data[[ii]]$variable_miss,
                                                        P = cb_obj$cb_model_para$P)
                      flag <- flag + 1
                    }
                    res <- Reduce(pmax, res)
                    specific_cs_purity <- lapply(1:3, function(ii){
                      specific_cs_purity[[ii]][i,j] <- specific_cs_purity[[ii]][j,i] <- res[ii]
                      return(specific_cs_purity[[ii]])
                    })
                  }
                }
                names(specific_cs_purity) <- c("min_abs_cor", "max_abs_cor", "median_abs_cor")
                
              } else {
                specific_cs_purity <- out_ucos$purity_each
                rownames(specific_cs_purity) <- specific_cs_names
              }
              
              # - cos&ucos purity
              cos <- cb_output$cos_details$cos$cos_index
              ncos <- length(cos)
              if (ncos!=0){
                empty_matrix <- matrix(NA, ncos, nucos)
                colnames(empty_matrix) <- specific_cs_names
                rownames(empty_matrix) <- names(cos)
                cos_ucos_purity <- lapply(1:3, function(ii) empty_matrix )
                for (i in 1:ncos){
                  for (j in 1:nucos){
                    cset1 <- cos[[i]]
                    cset2 <- specific_css$ucos_index[[j]]
                    y.i <- cb_output$cos_details$cos_outcomes$outcome_index[[i]]
                    y.j <- specific_outcomes$outcome_index[[j]]
                    yy <- unique(y.i, y.j)
                    res <- list()
                    flag <- 1
                    for (ii in yy){
                      X_dict <- cb_obj$cb_data$dict[ii]
                      res[[flag]] <- get_between_purity(cset1, cset2, X = cb_obj$cb_data$data[[X_dict]]$X,
                                                        Xcorr = cb_obj$cb_data$data[[X_dict]]$XtX,
                                                        N = cb_obj$cb_data$data[[ii]]$N,
                                                        miss_idx = cb_obj$cb_data$data[[ii]]$variable_miss,
                                                        P = cb_obj$cb_model_para$P)
                      flag <- flag + 1
                    }
                    res <- Reduce(pmax, res)
                    cos_ucos_purity <- lapply(1:3, function(ii){
                      cos_ucos_purity[[ii]][i,j] <- res[ii]
                      return(cos_ucos_purity[[ii]])
                    })
                  }
                }
                names(cos_ucos_purity) <- c("min_abs_cor", "max_abs_cor", "median_abs_cor")
              } else {
                cos_ucos_purity <- NULL
              }
              
              
              # - save coloc_results
              specific_results <- list("ucos" = specific_css,
                                       "ucos_outcomes" = specific_outcomes,
                                       "ucos_weight" = specific_w,
                                       "ucos_top_variables" = specific_cs_hits,
                                       "ucos_purity" = specific_cs_purity,
                                       "cos_ucos_purity" = cos_ucos_purity,
                                       "ucos_outcomes_delta" = cs_change)
            }
            
        } else {
            specific_results <- NULL
        }

        # - cb_model_para
        cb_model_para$N <- as.numeric(unlist(cb_model_para$N))
        cb_model_para$variables <- variables

        ll <- list("ucos_details" = specific_results,
                   "cb_model" = cb_model,
                   "cb_model_para" = cb_model_para)
    } else {
        # - cb_model_para
        cb_model_para$N <- as.numeric(unlist(cb_model_para$N))
        cb_model_para$variables <- variables
        ll <- list("ucos_detials" = NULL,
                   "cb_model" = cb_model,
                   "cb_model_para" = cb_model_para)
    }

    return(ll)

}

#' @noRd
#' @keywords cb_post_inference
get_summary_table_fm <- function(cb_output, outcome_names = NULL, gene_name = NULL){
  
  specific_cs <- cb_output$ucos_details
  if (length(specific_cs$ucos$ucos_index) != 0){
    
    cs_outcome <- cb_output$data_info$outcome_info$outcome_names
    if (!is.null(outcome_names)){ cs_outcome <- outcome_names }
    pip <- as.numeric(cb_output$pip)
    
    summary_table <- matrix(NA, nrow = length(specific_cs$ucos$ucos_index), ncol = 9)
    colnames(summary_table) <- c("outcomes", "ucos_id", "purity", 
                                 "top_variable", "top_variable_pip", "n_variables", "ucos_index",
                                 "ucos_variables", "ucos_variables_pip")
    summary_table <- as.data.frame(summary_table)
    summary_table[,1] <- cs_outcome[unlist(specific_cs$ucos_outcomes$outcome_index)]
    summary_table[,2] <- names(specific_cs$ucos$ucos_index)
    summary_table[,3] <- as.numeric(diag(as.matrix(specific_cs$ucos_purity$min_abs_cor)))
    summary_table[,4] <- unlist(sapply(specific_cs$ucos$ucos_variables, function(tmp) tmp[1]))    
    summary_table[,5] <- sapply(specific_cs$ucos$ucos_index, function(tmp) max(pip[tmp]))                                  
    summary_table[,6] <- as.numeric(sapply(specific_cs$ucos$ucos_index, length))
    summary_table[,7] <- unlist(sapply(specific_cs$ucos$ucos_index, function(tmp) paste0(tmp, collapse = "; ")))
    summary_table[,8] <- unlist(sapply(specific_cs$ucos$ucos_variables, function(tmp) paste0(tmp, collapse = "; ")))
    summary_table[,9] <- unlist(sapply(specific_cs$ucos$ucos_index, function(tmp) paste0(pip[tmp], collapse = "; ")))
    if (!is.null(gene_name)){ summary_table$gene_name <- gene_name }
    
  } else {
    summary_table <- NULL
  }
  return(summary_table)
  
}



cos_pvalue_filter <- function(cos_results, data_info = NULL, pv_cutoff = 1e-4){
  
  if (is.null(data_info))
    data_info <- get_data_info(cb_obj)
  
  n_cos <- length(cos_results$cos_results$cos$cos_index)
  filtered_cos <- list("cos_index" = NULL,
                       "cos_variables" = NULL)
  filtered_traits <- list("outcome_index" = NULL,
                          "outcome_name" = NULL)
  filtered_cos_vcp <- list()
  pp_remove <- c()
  for (i in 1:n_cos){
    cos_tmp <- cos_results$cos_results$cos$cos_index[[i]]
    cos_trait <- cos_results$cos_results$cos_outcomes$outcome_index[[i]]
    minPV <- sapply(cos_trait, function(tmp){
      z <- data_info$z[[tmp]]
      z <- z[cos_tmp]
      pv <- pchisq(z^2, 1, lower.tail = FALSE)
      min(pv)
    })
    pp <- which(minPV < pv_cutoff)
    if (length(pp) == 0 | length(pp) == 1) {
      pp_remove <- c(pp_remove, i)
      next
    }
    filtered_cos$cos_index <- c(filtered_cos$cos_index, 
                                cos_results$cos_results$cos$cos_index[i])
    filtered_cos$cos_variables <- c(filtered_cos$cos_variables, 
                                    cos_results$cos_results$cos$cos_variables[i])
    filtered_cos_vcp <- c(filtered_cos_vcp, cos_results$cos_results$cos_vcp[i])
    if (length(pp) == length(cos_trait)){
      filtered_traits$outcome_index <- c(filtered_traits$outcome_index, 
                                         cos_results$cos_results$cos_outcomes$outcome_index[i])
      filtered_traits$outcome_name <- c(filtered_traits$outcome_name, 
                                        cos_results$cos_results$cos_outcomes$outcome_name[i])
    } else {
      filtered_traits$outcome_index <- c(filtered_traits$outcome_index, 
                                         list(cos_results$cos_results$cos_outcomes$outcome_index[[i]][pp]))
      filtered_traits$outcome_name <- c(filtered_traits$outcome_name, 
                                        list(cos_results$cos_results$cos_outcomes$outcome_name[[i]][pp]))
    }
  }
  if (length(filtered_cos$cos_index)==0){
    cos_results <- list("cos_results" = NULL, "vcp" = NULL)
  } else {
    cos_results$cos_results$cos <- filtered_cos
    names(filtered_traits$outcome_index) <- names(cos_results$cos_results$cos$cos_index)
    names(filtered_traits$outcome_name) <- names(cos_results$cos_results$cos$cos_index)
    cos_results$cos_results$cos_outcomes <- filtered_traits
    cos_results$cos_results$cos_vcp <- filtered_cos_vcp
    coloc_hits <- coloc_hits_variablenames <- coloc_hits_names <- c()
    for (i in 1:length(filtered_cos_vcp)){
      inw <- filtered_cos_vcp[[i]]
      pp1 <- which(inw == max(inw))
      coloc_hits <- c(coloc_hits, pp1)
      coloc_hits_variablenames <- c(coloc_hits_variablenames, data_info$variables[pp1])
      if (length(pp1)==1){
        coloc_hits_names <- c(coloc_hits_names, names(filtered_cos_vcp)[i])
      } else {
        coloc_hits_names <- c(coloc_hits_names, paste0(names(filtered_cos_vcp)[i], ".", 1:length(pp1)))
      }
    }
    coloc_hits <- data.frame("top_index" = coloc_hits, "top_variables" = coloc_hits_variablenames) 
    rownames(coloc_hits) <- coloc_hits_names
    cos_results$cos_results$cos_top_variables <- coloc_hits
    
    # - change purity
    if (length(pp_remove)!=0){
      cos_results$cos_results$cos_purity$min_abs_cor <- as.matrix(cos_results$cos_results$cos_purity$min_abs_cor)[-pp_remove, -pp_remove, drop=FALSE]
      cos_results$cos_results$cos_purity$median_abs_cor <- as.matrix(cos_results$cos_results$cos_purity$median_abs_cor)[-pp_remove, -pp_remove, drop=FALSE]
      cos_results$cos_results$cos_purity$max_abs_cor <- as.matrix(cos_results$cos_results$cos_purity$max_abs_cor)[-pp_remove, -pp_remove, drop=FALSE]
      cos_results$cos_results$cos_PPC <- cos_results$cos_results$cos_PPC[-pp_remove]
    }
    
  }
  return(cos_results)
  
}


colocboost_get_npc_configs <- function(cb_results, npc_cutoff = 0.7, alpha = 1.5, coverage = 0.95){
  
  if (is.null(cb_results$cos_details)){
    warnings("No colocalization results in this region!")
    return(cb_results)
  }
  if (npc_cutoff == 0){
    warnings("No changes in the original results for the most likely colocalization with npc_cutoff=0")
    return(cb_results)
  }
  message(paste("Extracting colocalization results with npc_cutoff =", npc_cutoff, ".\n",
                "For each CoS, keep the outcomes configurations that the npc_outcome >", npc_cutoff, "."))
  
  remove_cos <- function(cb_results, remove_idx = NULL){
    if (length(remove_idx)==0){
      return(cb_results)
    }
    ncos <- length(cos_details$cos_top_variables)
    if (length(remove_idx) == ncos){
      cb_results$vcp <- NULL
      cb_results$cos_details <- NULL
      cb_results <- c(cb_results, list(vcp = NULL, cos_details = NULL))
      return(cb_results)
    }
    cos_details <- cb_results$cos_details
    cos_details$cos_top_variables <- cos_details$cos_top_variables[-remove_idx, , drop = FALSE]
    cos_details$cos$cos_index <- cos_details$cos$cos_index[-remove_idx]
    cos_details$cos$cos_variables <- cos_details$cos$cos_variables[-remove_idx]
    cos_details$cos_outcomes$outcome_index <- cos_details$cos_outcomes$outcome_index[-remove_idx]
    cos_details$cos_outcomes$outcome_name <- cos_details$cos_outcomes$outcome_name[-remove_idx]
    cos_details$cos_vcp <- cos_details$cos_vcp[-remove_idx]
    cos_details$cos_weights <- cos_details$cos_weights[-remove_idx]
    cos_details$cos_npc <- cos_details$cos_npc[-remove_idx]
    cos_details$cos_purity$min_abs_cor <- as.matrix(cos_details$cos_purity$min_abs_cor)[-remove_idx, -remove_idx, drop=FALSE]
    cos_details$cos_purity$median_abs_cor <- as.matrix(cos_details$cos_purity$median_abs_cor)[-remove_idx, -remove_idx, drop=FALSE]
    cos_details$cos_purity$max_abs_cor <- as.matrix(cos_details$cos_purity$max_abs_cor)[-remove_idx, -remove_idx, drop=FALSE]
    vcp <- as.vector(1 - apply(1 - do.call(cbind, cos_details$cos_vcp),1,prod))
    names(vcp) <- cb_results$data_info$variables
    cb_results$vcp <- vcp
    cb_results$cos_details <- cos_details
    return(cb_results)
  }
  
  cos_details <- cb_results$cos_details

  coloc_outcome_index <- coloc_outcome <- list()
  colocset_names <- c()
  for (i in 1:length(cos_details$cos$cos_index)){
    cos_npc_config <-  cos_details$cos_normalization_evidence[[i]]
    npc_outcome <- cos_npc_config$npc_outcome
    pos_pass <- which(npc_outcome>=npc_cutoff)
    if (length(pos_pass) == 0){
       coloc_outcome_index[[i]] <- 0
       coloc_outcome[[i]] <- 0
       colocset_names[i] <- paste0("remove", i)
    } else {
        coloc_outcome_index[[i]] <- sort(cos_npc_config$outcomes_index[pos_pass])
        coloc_outcome[[i]] <- rownames(cos_npc_config)[pos_pass]
        colocset_names[i] <- paste0("cos", i, ":", paste0(paste0("y", coloc_outcome_index[[i]]), collapse = "_"))
        if (grepl("merged", names(cos_details$cos$cos_index)[i])) {
          colocset_names[i] <- paste0(colocset_names[i], ":merged")
        }
    }
    
  }
  names(coloc_outcome) <- names(coloc_outcome_index) <- colocset_names
  cos_details$cos_outcomes <- list("outcome_index" = coloc_outcome_index, "outcome_name" = coloc_outcome)
  
  # - VCP
  cos_weights <- lapply(1:length(cos_details$cos_outcomes$outcome_index), function(idx){
    w <- cos_details$cos_weights[[idx]]
    config_idx <- cos_details$cos_outcomes$outcome_index[[idx]]
    if (length(config_idx)==1){
      if (config_idx == 0){
        return(matrix(1, nrow=nrow(cos_details$cos_weights[[idx]])))
      }
    }
    w_outcome <- colnames(w)
    config_outcome <- paste0("outcome", config_idx)
    pos <- which(w_outcome %in% config_outcome)
    w[, pos, drop=FALSE]
  })
  cos_details$cos_weights <- cos_weights
  int_weight <- lapply(cos_weights, get_integrated_weight, alpha = alpha)
  names(int_weight) <- names(cos_weights) <- colocset_names
  vcp <- as.vector(1 - apply(1 - do.call(cbind, int_weight),1,prod))
  names(vcp) <- cb_results$data_info$variables
  cb_results$vcp <- vcp
  cos_details$cos_vcp = int_weight
  
  # - resummary results
  cos_re_idx <- lapply(int_weight, function(w){ unlist(get_in_csets(w, coverage = coverage))})
  cos_re_var <- lapply(cos_re_idx, function(idx){ cb_results$data_info$variables[idx] })
  coloc_csets <- list("cos_index" = cos_re_idx, "cos_variables" = cos_re_var)
  cos_details$cos <- coloc_csets
  
  # - hits variables in each csets
  coloc_hits <- coloc_hits_variablenames <- coloc_hits_names <- c()
  for (i in 1:length(int_weight)){
    inw <- int_weight[[i]]
    if (length(unique(inw))==1){
      coloc_hits <- c(coloc_hits, 0)
      coloc_hits_names <- c(coloc_hits_names, names(int_weight)[i])
    } else {
      pp <- which(inw == max(inw))
      coloc_hits <- c(coloc_hits, pp)
      coloc_hits_variablenames <- c(coloc_hits_variablenames, cb_results$data_info$variables[pp])
      if (length(pp)==1){
        coloc_hits_names <- c(coloc_hits_names, names(int_weight)[i])
      } else {
        coloc_hits_names <- c(coloc_hits_names, paste0(names(int_weight)[i], ".", 1:length(pp)))
      }
    }
  }
  coloc_hits <- data.frame("top_index" = coloc_hits, "top_variables" = coloc_hits_variablenames) 
  rownames(coloc_hits) <- coloc_hits_names
  cos_details$cos_top_variables <- coloc_hits
  cb_results$cos_details <- cos_details
  
  
  # remove CoS does not include outcomes pass npc_cutoff
  remove <- grep("remove", colocset_names)
  cb_results <- remove_cos(cb_results, remove_idx = remove)
  cos_details <- cb_results$cos_details
  
  # remove CoS only with one trait
  n_outcome <- sapply(cos_details$cos_outcomes$outcome_index, length)
  single <- which(n_outcome == 1)
  if (length(single)==length(n_outcome)){
    
    # - all remaining the single outcome
    cb_results$cos_details <- cb_results$vcp <- NULL
    cb_results <- c(cb_results, list(vcp = NULL, cos_details = NULL))
    
  } else if (length(single)!=0 & length(single)!=length(n_outcome)){
    
    # - partial remaining the single outcome
    ucos_from_cos <- list("ucos" = list("ucos_index" = cos_details$cos$cos_index[single],
                                        "ucos_variables" = cos_details$cos$cos_variables[single]),
                          "ucos_outcomes" = list("outcome_idx" = cos_details$cos_outcomes$outcome_index[single],
                                                 "outcome_name" = cos_details$cos_outcomes$outcome_name[single]),
                          "ucos_weight" = cos_details$cos_weights[single],
                          "ucos_purity" = list("min_abs_cor" = as.matrix(cos_details$cos_purity$min_abs_cor)[single, single, drop=FALSE],
                                               "median_abs_cor" = as.matrix(cos_details$cos_purity$median_abs_cor)[single, single, drop=FALSE],
                                               "max_abs_cor" = as.matrix(cos_details$cos_purity$max_abs_cor)[single, single, drop=FALSE]),
                          "ucos_top_variables" = cos_details$cos_top_variables[single, , drop = FALSE])
    cb_results$ucos_from_cos <- ucos_from_cos
    cb_results <- remove_cos(cb_results, remove_idx = single)
    
  }
  
  # - refine and output
  target_outcome <- NULL
  target_idx <- which(cb_results$data_info$outcome_info$is_target)
  if (length(target_idx)!=0){
    target_outcome <- cb_results$data_info$outcome_info$outcome_names[target_idx]
  }
  cb_results$cos_summary <- get_cos_summary(cb_results, target_outcome = target_outcome)
  
  return(cb_results)
  
}














