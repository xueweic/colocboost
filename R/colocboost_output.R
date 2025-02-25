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
    variables <- cb_obj$cb_data$snp.names
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
get_cos_details <- function(cb_obj, coloc_out, data_info = NULL){

    if (is.null(data_info))
        data_info <- get_data_info(cb_obj)


    ### ----- Define the colocalization results
    coloc_sets <- coloc_out$csets
    if (length(coloc_sets)!=0){

        # - colocalized outcomes
        analysis_outcome <- cb_obj$cb_model_para$outcome_names
        coloc_outcome_index <- coloc_outcome <- list()
        colocset_names <- c()
        for (i in 1:length(coloc_out$csets)){
            cc <- coloc_out$avWeight[[i]]
            tmp_names <- colnames(cc)
            coloc_outcome_index[[i]] <- as.numeric(gsub("[^0-9.]+", "", tmp_names))
            coloc_outcome[[i]] <- analysis_outcome[coloc_outcome_index[[i]]]
            colocset_names[i] <- paste0("cos", i, ":", paste0(paste0("y", coloc_outcome_index[[i]]), collapse = "_"))
            if (grepl("Merge", names(coloc_sets)[i])) {
              colocset_names[i] <- paste0(colocset_names[i], ":merge")
            }
        }
        names(coloc_outcome) <- names(coloc_outcome_index) <- colocset_names
        coloc_outcomes <- list("outcome_index" = coloc_outcome_index, "outcome_name" = coloc_outcome)

        # - colocalized sets for variables
        coloc_csets_snpidx <- coloc_out$csets
        coloc_csets_snpnames <- lapply(coloc_csets_snpidx, function(coloc_tmp){ cb_obj$cb_model_para$variables[coloc_tmp] })
        coloc_csets_snpidx <- lapply(coloc_csets_snpnames, function(snp) match(snp, data_info$variables))
        names(coloc_csets_snpidx) <- names(coloc_csets_snpnames) <- colocset_names
        coloc_csets_original <- list("cos_index" = coloc_csets_snpidx, "cos_variables" = coloc_csets_snpnames)

        # - colocalized set cs_change
        cs_change <- coloc_out$cs_change
        rownames(cs_change) <- colocset_names
        colnames(cs_change) <- analysis_outcome

        # - VCP
        int_weight <- lapply(coloc_out$avWeight, get_integrated_weight, alpha = cb_obj$cb_model_para$alpha)
        names(int_weight) <- colocset_names
        int_weight <- lapply(int_weight, function(inw) {
            pos <- match(data_info$variables, cb_obj$cb_model_para$variables)
            return(inw[pos])
        })
        vcp <- as.vector(1 - apply(1 - do.call(cbind, int_weight),1,prod))
        names(vcp) <- data_info$variables
        
        # - resummary results
        cos_re_idx <- lapply(int_weight, function(w){ unlist(get_in_csets(w, coverage = cb_obj$cb_model_para$coverage))})
        cos_re_var <- lapply(cos_re_idx, function(idx){ data_info$variables[idx] })
        coloc_csets <- list("cos_index" = cos_re_idx, "cos_variables" = cos_re_var)
        
        # - hits variables in each csets
        coloc_hits <- coloc_hits_snpnames <- coloc_hits_names <- c()
        for (i in 1:length(int_weight)){
            inw <- int_weight[[i]]
            pp <- which(inw == max(inw))
            coloc_hits <- c(coloc_hits, pp)
            coloc_hits_snpnames <- c(coloc_hits_snpnames, data_info$variables[pp])
            if (length(pp)==1){
                coloc_hits_names <- c(coloc_hits_names, names(int_weight)[i])
            } else {
                coloc_hits_names <- c(coloc_hits_names, paste0(names(int_weight)[i], ".", 1:length(pp)))
            }
        }
        coloc_hits <- data.frame("top_index" = coloc_hits, "top_variables" = coloc_hits_snpnames) 
        rownames(coloc_hits) <- coloc_hits_names

        # - purity
        ncos <- length(coloc_sets)
        if (ncos >= 2){
            empty_matrix <- matrix(NA, ncos, ncos)
            colnames(empty_matrix) <- rownames(empty_matrix) <- colocset_names
            csets_purity <- lapply(1:3, function(ii){
              diag(empty_matrix) <- coloc_out$purity[,ii]
              return(empty_matrix)
            })
            for (i in 1:(ncos-1)){
                for (j in (i+1):ncos){
                    cset1 <- coloc_sets[[i]]
                    cset2 <- coloc_sets[[j]]
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
                                                          miss_idx = cb_obj$cb_data$data[[ii]]$snp_miss,
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
                              "cos_purity" = csets_purity,
                              "cos_top_variables" = coloc_hits,
                              "cos_outcomes_delta" = cs_change)
        
        
        # - missing snp and warning message
        missing_snps_idx <- Reduce(union, lapply(cb_obj$cb_data$data,function(cb) cb$snp_miss ))
        missing_snps <- cb_obj$cb_model_para$variables[missing_snps_idx]
        cos_missing_snps_idx <-  lapply(coloc_csets_original$cos_variables, function(snp){
          missing <- intersect(snp, missing_snps)
          if (length(missing)!=0){
            match(missing, data_info$variables)
          } else { NULL }
        })
        cos_missing_snps <- lapply(cos_missing_snps_idx, function(snp){if(!is.null(snp)){data_info$variables[snp]} else {NULL}})
        warning_needed <- any(!sapply(cos_missing_snps, is.null))
        if (warning_needed){
           is_missing <- which(!sapply(cos_missing_snps, is.null))
           cos_missing_snps_idx <- cos_missing_snps_idx[is_missing]
           cos_missing_snps <- cos_missing_snps[is_missing]
           cos_missing_vcp <- lapply(cos_missing_snps_idx, function(idx){ vcp[idx] })
           warning_message <- paste("CoS", paste(names(cos_missing_snps_idx), collapse = ","),
                                    "contains missing variables in at least one outcome.",
                                    "The missing variables will cause the ~0 VCP scores.")
           cos_warnings <- list("cos_missing_info" = list("index" = cos_missing_snps_idx,
                                                          "variables" = cos_missing_snps,
                                                          "vcp" = cos_missing_vcp),
                                "warning_message" = warning_message)
           coloc_results$cos_warnings <- cos_warnings
        }

    } else {
        coloc_results <- NULL
        vcp <- NULL
    }
    return(list("coloc_results"=coloc_results, "vcp"=vcp))

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
get_cos_summary <- function(cb_output, outcome_names = NULL, gene_name = NULL, target_outcome = NULL){

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
    
        summary_table <- matrix(NA, nrow = length(coloc_sets), ncol = 10)
        colnames(summary_table) <- c("target_outcome", "colocalized_outcomes", "cos_id", "purity", 
                                     "top_variable", "top_variable_vcp", "n_variables", "colocalized_index",
                                     "colocalized_variables", "colocalized_variables_vcp")
        summary_table <- as.data.frame(summary_table)
        summary_table[,1] <- FALSE
        summary_table[,2] <- unlist(sapply(coloc_outcome, function(tmp) paste0(tmp, collapse = "; ")))
        summary_table[,3] <- names(coloc_sets)
        summary_table[,4] <- as.numeric(diag(as.matrix(cb_output$cos_details$cos_purity$min_abs_cor)))
        summary_table[,5] <- unlist(sapply(cb_output$cos_details$cos$cos_variables, function(tmp) tmp[1]))    
        summary_table[,6] <- sapply(coloc_sets, function(tmp) max(vcp[tmp]))                                  
        summary_table[,7] <- as.numeric(sapply(coloc_sets, length))
        summary_table[,8] <- unlist(sapply(coloc_sets, function(tmp) paste0(tmp, collapse = "; ")))
        summary_table[,9] <- unlist(sapply(cb_output$cos_details$cos$cos_variables, function(tmp) paste0(tmp, collapse = "; ")))
        summary_table[,10] <- unlist(sapply(coloc_sets, function(tmp) paste0(vcp[tmp], collapse = "; ")))
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
get_full_output <- function(cb_obj, past_out = NULL, variables = NULL, cb_output = NULL){

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
        out_single <- past_out$single
        # - single sets
        if (!is.null(out_single$csets_each)){
            out_single$csets_each <- lapply(out_single$csets_each, function(cs){
                match(cb_model_para$variables[cs], variables)
            })
            out_single$avW_csets_each <- out_single$avW_csets_each[ordered,,drop=FALSE]
            
            # - re-orginize specific results
            analysis_outcome <- cb_obj$cb_model_para$outcome_names
            specific_outcome_index <- specific_outcome <- list()
            specific_cs_names <- c()
            for (i in 1:length(out_single$csets_each)){
              cc <- out_single$avW_csets_each[,i,drop=FALSE] 
              tmp_names <- colnames(cc)
              specific_outcome_index[[i]] <- as.numeric(gsub(".*Y([0-9]+).*", "\\1", tmp_names))
              specific_outcome[[i]] <- analysis_outcome[specific_outcome_index[[i]]]
              specific_cs_names[i] <- paste0("ucos", i, ":y", specific_outcome_index[[i]])
            }
            names(specific_outcome) <- names(specific_outcome_index) <- specific_cs_names
            specific_outcomes <- list("outcome_index" = specific_outcome_index, "outcome_name" = specific_outcome)
            
            # - specific sets for variables
            specific_cs_snpidx <- out_single$csets_each
            specific_cs_snpnames <- lapply(specific_cs_snpidx, function(specific_tmp){ cb_obj$cb_model_para$variables[specific_tmp] })
            specific_cs_snpidx <- lapply(specific_cs_snpnames, function(snp) match(snp, variables))
            names(specific_cs_snpidx) <- names(specific_cs_snpnames) <- specific_cs_names
            specific_css <- list("ucos_index" = specific_cs_snpidx, "ucos_variables" = specific_cs_snpnames)
            
            # - specific set cs_change
            cs_change <- out_single$change_obj_each
            rownames(cs_change) <- specific_cs_names
            colnames(cs_change) <- analysis_outcome
            index_change <- as.data.frame(which(cs_change != 0, arr.ind = TRUE))
            change_outcomes <- analysis_outcome[index_change$col]
            change_values <- diag(as.matrix(cs_change[index_change$row, index_change$col]))
            cs_change <- data.frame("ucos_outcome" = change_outcomes, "ucos_delta" = change_values)
            
            # - ucos_weight
            specific_w <- lapply(1:ncol(out_single$avW_csets_each), function(ii) out_single$avW_csets_each[,ii,drop=FALSE])
            names(specific_w) <- specific_cs_names
            
            # - hits variables in each csets
            cs_hits <- sapply(1:length(specific_w), function(jj){ inw=specific_w[[jj]]; sample(which(inw == max(inw)),1) })
            cs_hits_snpnames <- sapply(cs_hits, function(ch) variables[ch] )
            specific_cs_hits <- data.frame("top_index" = cs_hits, "top_variables" = cs_hits_snpnames)   # save
            rownames(specific_cs_hits) <- specific_cs_names
            
            # - purity
            nucos <- length(out_single$csets_each)
            if (nucos >= 2){
              empty_matrix <- matrix(NA, nucos, nucos)
              colnames(empty_matrix) <- rownames(empty_matrix) <- specific_cs_names
              specific_cs_purity <- lapply(1:3, function(ii){
                diag(empty_matrix) <- out_single$purity_each[,ii]
                return(empty_matrix)
              })
              for (i in 1:(nucos-1)){
                for (j in (i+1):nucos){
                  cset1 <- out_single$csets_each[[i]]
                  cset2 <- out_single$csets_each[[j]]
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
                                                      miss_idx = cb_obj$cb_data$data[[ii]]$snp_miss,
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
              specific_cs_purity <- out_single$purity_each
              rownames(specific_cs_purity) <- specific_cs_names
            }
            
            # - cos&ucos purity
            cos <- cb_output$cos_details$cos$cos_index
            ucos <- 
            ncos <- length(cos)
            if (ncos!=0){
              empty_matrix <- matrix(NA, ncos, nucos)
              colnames(empty_matrix) <- specific_cs_names
              rownames(empty_matrix) <- names(cos)
              cos_ucos_purity <- lapply(1:3, function(ii) empty_matrix )
              for (i in 1:ncos){
                for (j in 1:nucos){
                  cset1 <- cos[[i]]
                  cset2 <- specific_cs_snpidx[[j]]
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
                                                      miss_idx = cb_obj$cb_data$data[[ii]]$snp_miss,
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
    summary_table[,3] <- as.numeric(diag(specific_cs$ucos_purity$min_abs_cor))
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

