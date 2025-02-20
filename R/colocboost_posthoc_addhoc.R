

merge_coloc_single <- function(cb_obj, out_coloc, out_single, coverage = 0.95,
                               min_abs_corr = 0.5, tol = 1e-9,
                               between_purity = 0.8){


    change_obj_each <- out_single$change_obj_each
    coloc_sets <- out_coloc$csets$csets
    csets_each <- out_single$csets_each

    # - remove overlap between coloc_sets and single_sets
    is_overlap <- is_highLD <- c()
    for (i in 1:length(coloc_sets)){
        # cset1 <- coloc_sets[[i]]
        coloc_trait <- out_coloc$csets$coloc_traits[[i]]
        ttmp <- as.numeric(gsub("[^0-9.]+", "", colnames(out_coloc$csets$avWeight[[i]])))
        pos_name <- order(ttmp)
        out_coloc$csets$avWeight[[i]] <- out_coloc$csets$avWeight[[i]][,pos_name]
        for (j in 1:length(csets_each)){
            cset2 <- csets_each[[j]]
            fine_trait <- out_single$single_trait[j]

            if (fine_trait %in% coloc_trait){

                # - if fine_y in coloc_y, we check only the overlap
                pp <- which(coloc_trait == fine_trait)
                w_tmp <- out_coloc$csets$avWeight[[i]][,pp]
                cset1 <- get_in_csets(w_tmp, coverage = coverage)[[1]]
                # - addhoc: if between_purity > 0.8, remove single sets
                X_dict <- cb_obj$cb_data$dict[fine_trait]
                res <- get_between_purity(cset1, cset2, X = cb_obj$cb_data$data[[X_dict]]$X,
                                          Xcorr = cb_obj$cb_data$data[[X_dict]]$XtX,
                                          N = cb_obj$cb_data$data[[fine_trait]]$N,
                                          miss_idx = cb_obj$cb_data$data[[fine_trait]]$snp_miss,
                                          P = cb_obj$cb_model_para$P)
                # is.between <- length(intersect(cset1, cset2)) != 0
                is.between <- (abs(res[2]-1) < tol)
                if (is.between){
                    is_overlap = c(is_overlap, j)
                } else {
                    is.higLD <- res[2] > between_purity
                    if (is.higLD){ is_highLD = c(is_highLD, j)}
                }

            } else if (!(fine_trait %in% coloc_trait)){

                # - if fine_y not in coloc_y, we check overlap and also min_purity
                change_obj_coloc <- out_coloc$csets$cs_change
                cset1 <- coloc_sets[[i]]
                res <- list()
                for (ii in 1:cb_obj$cb_model_para$L){
                    X_dict <- cb_obj$cb_data$dict[ii]
                    res[[ii]] <- get_between_purity(cset1, cset2, X = cb_obj$cb_data$data[[X_dict]]$X,
                                                    Xcorr = cb_obj$cb_data$data[[X_dict]]$XtX,
                                                    N = cb_obj$cb_data$data[[ii]]$N,
                                                    miss_idx = cb_obj$cb_data$data[[ii]]$snp_miss,
                                                    P = cb_obj$cb_model_para$P)
                }
                res <- Reduce(pmax, res)
                min_between <- res[1]
                max_between <- res[2]
                ave_between <- res[3]
                # is.between <- ((min_between>min_abs_corr) & (abs(max_between-1) < tol))
                is.between <- (min_between>min_abs_corr) * (abs(max_between-1)<tol) * (ave_between>between_purity)
                if (is.between){
                    is_overlap = c(is_overlap, j)
                    # -- add weight
                    add_avW <- out_single$avW_csets_each[,j]
                    out_coloc$csets$avWeight[[i]] <- cbind(add_avW, out_coloc$csets$avWeight[[i]])
                    colnames(out_coloc$csets$avWeight[[i]])[1] <- paste0("trait", fine_trait)
                    ttmp <- as.numeric(gsub("[^0-9.]+", "", colnames(out_coloc$csets$avWeight[[i]])))
                    pos_name <- order(ttmp)
                    out_coloc$csets$avWeight[[i]] <- out_coloc$csets$avWeight[[i]][,pos_name]
                    change_obj_coloc_i <- change_obj_coloc[i,]
                    change_obj_each_j <- change_obj_each[j,]
                    out_coloc$csets$cs_change[i,] <- pmax(change_obj_coloc_i, change_obj_each_j)
                    coloc_trait <- c(coloc_trait, fine_trait)
                    out_coloc$csets$coloc_traits[[i]] <- coloc_trait
                }

            }
        }
    }

    is_overlap <- unique(c(is_overlap,is_highLD))
    if (length(is_overlap) != 0){
        if (length(is_overlap) == length(csets_each)){
            out_single$csets_each = NULL
            out_single$avW_csets_each <- NULL
            out_single$change_obj_each <- NULL
            out_single$purity_each <- NULL
        } else {
            out_single$csets_each = csets_each[-is_overlap]
            out_single$avW_csets_each <- out_single$avW_csets_each[,-is_overlap,drop = FALSE]
            out_single$change_obj_each <- change_obj_each[-is_overlap,,drop = FALSE]
            out_single$purity_each <- out_single$purity_each[-is_overlap,,drop = FALSE]
        }
    }
    if (length(out_single$csets_each) == 0){
        out_single$csets_each = NULL
        out_single$avW_csets_each <- NULL
        out_single$change_obj_each <- NULL
        out_single$purity_each <- NULL
    }
    ll <- list("single" = out_single, "coloc" = out_coloc)
    return(ll)

}




merge_coloc_single_backup <- function(cb_obj, out_coloc, out_single,
                                      min_abs_corr = 0.5,
                                      tol = 1e-9){


    change_obj_coloc <- out_coloc$csets$cs_change
    change_obj_each <- out_single$change_obj_each
    coloc_sets <- out_coloc$csets$csets
    csets_each <- out_single$csets_each

    # - remove overlap between coloc_sets and single_sets
    is_overlap <- c()
    for (i in 1:length(coloc_sets)){
        cset1 <- coloc_sets[[i]]
        coloc_trait <- out_coloc$csets$coloc_traits[[i]]
        ttmp <- as.numeric(gsub("[^0-9.]+", "", colnames(out_coloc$csets$avWeight[[i]])))
        pos_name <- order(ttmp)
        out_coloc$csets$avWeight[[i]] <- out_coloc$csets$avWeight[[i]][,pos_name]
        for (j in 1:length(csets_each)){
            cset2 <- csets_each[[j]]
            fine_trait <- out_single$single_trait[j]
            res <- list()
            for (ii in 1:cb_obj$cb_model_para$L){
                X_dict <- cb_obj$cb_data$dict[ii]
                res[[ii]] <- get_between_purity(cset1, cset2, X = cb_obj$cb_data$data[[X_dict]]$X,
                                                Xcorr = cb_obj$cb_data$data[[X_dict]]$XtX,
                                                N = cb_obj$cb_data$data[[ii]]$N,
                                                miss_idx = cb_obj$cb_data$data[[ii]]$snp_miss,
                                                P = cb_obj$cb_model_para$P)
            }
            res <- Reduce(pmax, res)
            min_between <- res[1]
            max_between <- res[2]
            # is.between <- (min_between>min_abs_corr) * (abs(max_between-1) < tol)
            is.between <-  abs(max_between-1) < tol
            if (is.between){
                # if (all(cset1 %in% cset2) | all(cset2 %in% cset1)){
                is_overlap = c(is_overlap, j)
                # - if trait difference, we need to merge it
                if ( !(fine_trait %in% coloc_trait) ){
                    # -- add weight
                    add_avW <- out_single$avW_csets_each[,j]
                    out_coloc$csets$avWeight[[i]] <- cbind(add_avW, out_coloc$csets$avWeight[[i]])
                    colnames(out_coloc$csets$avWeight[[i]])[1] <- paste0("trait", fine_trait)
                    ttmp <- as.numeric(gsub("[^0-9.]+", "", colnames(out_coloc$csets$avWeight[[i]])))
                    pos_name <- order(ttmp)
                    out_coloc$csets$avWeight[[i]] <- out_coloc$csets$avWeight[[i]][,pos_name]
                    change_obj_coloc_i <- change_obj_coloc[i,]
                    change_obj_each_j <- change_obj_each[j,]
                    out_coloc$csets$cs_change[i,] <- pmax(change_obj_coloc_i, change_obj_each_j)
                    out_coloc$csets$coloc_traits[[i]] <- c(coloc_trait, fine_trait)
                }
            }
        }
    }

    is_overlap <- unique(is_overlap)
    if (length(is_overlap) != 0){
        if (length(is_overlap) == length(csets_each)){
            out_single$csets_each = NULL
            out_single$avW_csets_each <- NULL
            out_single$change_obj_each <- NULL
            out_single$purity_each <- NULL
        } else {
            out_single$csets_each = csets_each[-is_overlap]
            out_single$avW_csets_each <- out_single$avW_csets_each[,-is_overlap,drop = FALSE]
            out_single$change_obj_each <- change_obj_each[-is_overlap,]
            out_single$purity_each <- out_single$purity_each[-is.overlap,]
        }
    }
    if (length(out_single$csets_each) == 0){
        out_single$csets_each = NULL
        out_single$avW_csets_each <- NULL
        out_single$change_obj_each <- NULL
        out_single$purity_each <- NULL
    }
    ll <- list("single" = out_single, "coloc" = out_coloc)
    return(ll)

}




merge_single <- function(cb_obj, past_out,
                         min_abs_corr = 0.5,
                         median_abs_corr = NULL,
                         n_purity = 100,
                         between_purity = 0.8,
                         tol = 1e-9){


    out_single <- past_out$single
    out_coloc <- past_out$coloc
    csets_each <- out_single$csets_each
    change_obj_each <- out_single$change_obj_each

    # calculate between purity
    ncsets <- length(csets_each)
    min_between <- max_between <- ave_between <- matrix(0, nrow =  ncsets, ncol =  ncsets)
    for (i.between in 1:(ncsets-1)){
        for (j.between in (i.between+1):ncsets){
            cset1 <- csets_each[[i.between]]
            cset2 <- csets_each[[j.between]]
            y.i <- out_single$single_trait[i.between]
            y.j <- out_single$single_trait[j.between]
            if (y.i == y.j){ next }
            yy <- c(y.i, y.j)
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
            if (min_abs_corr == 0){
              res <- Reduce(pmin, res)
            } else {
              res <- Reduce(pmax, res)
            }
            min_between[i.between, j.between] <- min_between[j.between, i.between] <- res[1]
            max_between[i.between, j.between] <- max_between[j.between, i.between] <- res[2]
            ave_between[i.between, j.between] <- ave_between[j.between, i.between] <- res[3]
        }
    }
    is.between <- (min_between>min_abs_corr) * (abs(max_between-1)<tol) * (ave_between>between_purity)
    if (sum(is.between) != 0){
        temp <- sapply(1:nrow(is.between), function(x){
            tt <- c(x, which(is.between[x,] != 0))
            return(paste0(sort(tt), collapse = ";"))
        })
        temp <- merge_sets(temp)
        potential_merged <- lapply(temp, function(x) as.numeric(unlist(strsplit(x, ";"))))
        potential_merged <- potential_merged[which(sapply(potential_merged, length) >= 2)]
        coloc_sets_merged <- avWeight_merged <-
            cs_change_merged <- coloc_traits_merged <- list()
        is_merged <- c()
        for (i.m in 1:length(potential_merged)){
            temp_set <- as.numeric(potential_merged[[i.m]])
            is_merged <- c(is_merged, temp_set)
            # define merged set
            coloc_sets_merged <- c(coloc_sets_merged, list( unique(unlist(csets_each[temp_set])) ))
            # refine avWeight
            merged <- out_single$avW_csets_each[, temp_set]
            unique_coloc_traits <- as.numeric(gsub(".*Y(\\d+).*", "\\1", colnames(merged)))
            colnames(merged) <- paste0("trait", unique_coloc_traits)
            coloc_traits_merged <- c(coloc_traits_merged,
                                     list(unique_coloc_traits[order(unique_coloc_traits)]))
            # colnames(temp) <- unique_coloc_traits
            avWeight_merged <- c(avWeight_merged, list(merged))
            # refine cs_change
            change_cs_tmp <- change_obj_each[temp_set, , drop = FALSE]
            cs_change_merged <- c(cs_change_merged,
                                  list(apply(change_cs_tmp, 2, max)))
        }

        if (length(is_merged) != 0){

            # --- check merged coloc set purity
            purity = vector(mode='list', length=length(coloc_sets_merged))
            for(ee in 1:length(coloc_sets_merged)){
                coloc_t <- coloc_traits_merged[[ee]]
                p_tmp <- c()
                for (i3 in coloc_t){
                    pos <- coloc_sets_merged[[ee]]
                    X_dict <- cb_obj$cb_data$dict[i3]
                    if (!is.null(cb_obj$cb_data$data[[X_dict]]$XtX)){
                        pos <- match(pos, setdiff(1:cb_obj$cb_model_para$P, cb_obj$cb_data$data[[i3]]$snp_miss))
                        # - it could happen since merging
                        if(sum(is.na(pos)) != 0){pos <- as.numeric(na.omit(pos))}
                    }
                    tmp <- matrix(get_purity(pos,X=cb_obj$cb_data$data[[X_dict]]$X,
                                             Xcorr=cb_obj$cb_data$data[[X_dict]]$XtX,
                                             N=cb_obj$cb_data$data[[i3]]$N,n=n_purity),1,3)
                    p_tmp <- rbind(p_tmp, tmp)
                }
                purity[[ee]] <- matrix(apply(p_tmp, 2, max),1,3)

            }
            purity_all <- do.call(rbind, purity)
            purity_all = as.data.frame(purity_all)
            colnames(purity_all) = c("min_abs_corr","mean_abs_corr","median_abs_corr")
            if (is.null(median_abs_corr)) {
                is_pure = which(purity_all[,1] >= min_abs_corr)
            } else {
                is_pure = which(purity_all[,1] >= min_abs_corr | purity_all[,3] >= median_abs_corr)
            }
            if (length(is_pure) > 0){
                # add coloc
                merged_colocsets <- coloc_sets_merged[is_pure]
                names(merged_colocsets) <- paste0("MergeCS", 1:length(is_pure))
                out_coloc$csets$csets <- c(out_coloc$csets$csets,
                                           merged_colocsets)
                out_coloc$csets$purity <- rbind(out_coloc$csets$purity,
                                                purity_all[is_pure, ])
                out_coloc$csets$evidence_strength <- c(out_coloc$csets$evidence_strength,
                                                       rep(0.95,length(is_pure)))
                cs_change_tmp <- do.call(rbind, cs_change_merged)[is_pure, , drop=FALSE]
                colnames(cs_change_tmp) <- paste0("change_obj_", 1:cb_obj$cb_model_para$L)
                rownames(cs_change_tmp) <- names(merged_colocsets)
                out_coloc$csets$cs_change <- rbind(out_coloc$csets$cs_change,
                                                   cs_change_tmp)
                out_coloc$csets$avWeight <- c(out_coloc$csets$avWeight,
                                              avWeight_merged[is_pure])
                out_coloc$csets$coloc_traits <- c(out_coloc$csets$coloc_traits,
                                                      coloc_traits_merged[is_pure])
            }

            # - remove single
            if (length(is_merged) == length(csets_each)){
                out_single$csets_each <- NULL
                out_single$avW_csets_each <- NULL
                out_single$change_obj_each <- NULL
                out_single$purity_each <- NULL
            } else {
                out_single$csets_each <- csets_each[-is_merged]
                out_single$avW_csets_each <- out_single$avW_csets_each[, -is_merged, drop = FALSE]
                out_single$change_obj_each <- change_obj_each[-is_merged, , drop = FALSE]
                out_single$purity_each <- out_single$purity_each[-is_merged, , drop = FALSE]
            }
        }
    }
    ll <- list("single" = out_single, "coloc" = out_coloc)
    return(ll)
}


get_vcp <- function(past_out, P){
    if (!is.null(past_out$coloc$csets$csets)){
        avW_coloc_vcp <- sapply(past_out$coloc$csets$avWeight, get_integrated_weight)
    } else {avW_coloc_vcp <- NULL}
    # avW_csets_each <- past_out$single$avW_csets_each
    # all_weight <- cbind(avW_coloc_vcp, avW_csets_each) # out_p$csets$avWeight)
    all_weight <- avW_coloc_vcp
    if (length(all_weight) == P){
        all_weight <- as.matrix(unlist(all_weight))
    }
    if (!is.null(all_weight)){
        all_weight <- apply(all_weight, 2, as.numeric)
        all_weight <- as.matrix(all_weight)
        vcp <- as.vector(1 - apply(1 - all_weight,1,prod))
    } else {
        vcp <- rep(0, P)
    }
    return(vcp)
}


get_pip <- function(past_out, R, P){
    if (length(past_out$coloc$csets$csets)!=0){
        av_coloc <- do.call(cbind, past_out$coloc$csets$avWeight)
    } else {av_coloc = NULL}
    if (length(past_out$single$csets_each)!=0){
        av_noncoloc <- past_out$single$avW_csets_each
        tmp <- do.call(rbind, strsplit(colnames(av_noncoloc), ":"))
        colnames(av_noncoloc) <- paste0("trait", gsub("[^0-9.]+", "", tmp[,2]))
    } else {av_noncoloc = NULL}
    av_all <- cbind(av_coloc, av_noncoloc)
    pip <- vector(mode='list', length=R)
    if (!is.null(av_all)){
        av_name <- colnames(av_all)
        for (i in 1:R){
            pos <- grep(i, av_name)
            if (length(pos) != 0){
                av_i <- as.matrix(av_all[, pos])
                pip[[i]] <- as.vector(1 - apply(1-av_i,1,prod))
            } else {
                pip[[i]] <- rep(0,P)
            }
        }
    }
    return(pip)
}

