

get_integrated_weight <- function(avWeight, func_intw = "fun_R", alpha = 1.5){
    
    if (func_intw == "prod"){
        av <- apply(avWeight, 1, function(w) prod(w))
    } else if (func_intw == "sqrt_prod"){
        av <- apply(avWeight, 1, function(w) prod(sqrt(w)))
    } else if (func_intw == "fun_R"){
        av <- apply(avWeight, 1, function(w) prod(w^(alpha/ncol(avWeight))))
    } else {
        stop("Please check func_intw! ")
    }
    return(av / sum(av))
}


get_in_csets <- function(weights, coverage = 0.95){
    
    temp <- order(weights, decreasing=T)
    csets <- temp[1:min(which(cumsum(weights[temp]) > coverage))] # 95%
    return(list(csets))
    
}

# - Fast calculate correlation matrix
get_cormat <- function(X, intercepte = FALSE){
  X = t(X)
  # Center each variable
  if (!intercepte){
    X = X - rowMeans(X)
  }
  # Standardize each variable
  X = X / sqrt(rowSums(X^2))
  # Calculate correlations
  cr = tcrossprod(X)
  return(cr)
}

check_null_post <- function(cb_obj, 
                            coloc_sets_temp,
                            coloc_outcomes,
                            check_null = 0.1,
                            check_null_method = "profile",
                            weaker_ucos = TRUE){
    
    extract_last <- function(lst) { tail(lst, n = 1)}
    extract_first <- function(lst) { head(lst, n = 1)}
    
    get_profile <- function(cs_beta, X = NULL, Y = NULL, N = NULL,
                            XtX = NULL, YtY = NULL, XtY = NULL, miss_idx){
        
        if (!is.null(X)){
            mean((Y-X%*%as.matrix(cs_beta))^2) *N/(N-1)
        } else if (!is.null(XtY)){
            scaling_factor <- if (!is.null(N)) (N - 1) else 1
            yty <- YtY / scaling_factor
            xtx <- XtX
            if (length(miss_idx)!=0){
                xty <- XtY[-miss_idx] / scaling_factor
                cs_beta <- cs_beta[-miss_idx]
            } else { xty <- XtY  / scaling_factor  }
            
            yty - 2*sum(cs_beta*xty) + sum( (xtx %*% as.matrix(cs_beta)) * cs_beta )
        }
        
    }
    
    get_cs_obj <- function(cs_beta, res, tau, func_prior, lambda, adj_dep, LD_obj,
                           X = NULL, Y = NULL, N = NULL,
                           XtX = NULL, YtY = NULL, XtY = NULL, miss_idx = NULL){
        
        correlation <- get_correlation(X = X, res = res, XtY = XtY, N = N, YtY = YtY, 
                                       XtX = XtX, beta_k = cs_beta, miss_idx = miss_idx)
        z <- get_z(correlation, n=N, res)
        abs_cor <- abs(correlation)
        jk <- which(abs_cor == max(abs_cor))
        jk <- ifelse(length(jk) == 1, jk, sample(jk,1))
        P <- length(z)
        ld_jk <- get_LD_jk(jk, X = X, XtX = XtX, N = N,
                             remain_idx = setdiff(1:P, miss_idx), P = P)
        ld_feature <- sqrt(abs(ld_jk))
        # - calculate delta
        delta <- boost_KL_delta(z = z, ld_feature = ld_feature, adj_dep=adj_dep,
                                func_prior = func_prior, lambda = lambda)
        scaling_factor <- if (!is.null(N)) (N - 1) else 1
        cov_Xtr <- if(!is.null(X)){
            t(X)%*%res / scaling_factor
        } else {
            res / scaling_factor
        }
        obj_ld <- if (LD_obj) ld_feature else rep(1, length(ld_feature))
        obj_ld[miss_idx] <- 0
        exp_term <- adj_dep * obj_ld * (abs(cov_Xtr))
        return(tau * matrixStats::logSumExp(exp_term / tau + log(delta)))
        
    }
    
    update_res <- function(X = NULL, Y = NULL, XtX = NULL, XtY = NULL, N = NULL, cs_beta, miss_idx){
        if (!is.null(X)){
            return(Y - X%*%cs_beta)
        } else if (!is.null(XtX)){
            scaling.factor <- if (!is.null(N)) N-1 else 1
            xtx <- XtX / scaling.factor
            if (length(miss_idx)!=0){
                xty <- XtY[-miss_idx] / scaling.factor
                res.tmp <- rep(0, length(XtY))
                res.tmp[-miss_idx] <- xty - xtx%*%cs_beta[-miss_idx]
            } else {
                xty <- XtY / scaling.factor
                res.tmp <- xty - xtx%*%cs_beta
            }
            return(res.tmp)
        }
    }
    # - add hoc
    cut <- if(length(cb_obj$cb_data)==1) 0.2 else 1
    
    # ----- null filtering
    cb_data <- cb_obj$cb_data
    cs_change <- check_cs_change <- matrix(0, nrow = length(coloc_sets_temp), ncol = cb_obj$cb_model_para$L)
    colnames(cs_change) <- colnames(check_cs_change) <- paste0("change_obj_", 1:cb_obj$cb_model_para$L)
    max_change <- matrix(0, nrow = length(coloc_sets_temp), ncol = cb_obj$cb_model_para$L)
    for (i in 1:length(coloc_sets_temp)){
        cs_variants <- as.numeric(unlist(coloc_sets_temp[[i]]))
        for (j in coloc_outcomes){
            cs_beta <- cb_obj$cb_model[[j]]$beta
            cs_beta[cs_variants] <- 0
            X_dict <- cb_data$dict[j]
            adj_dep <- cb_data$data[[j]]$dependency
            if (check_null_method == "profile"){
                cs_profile <- get_profile(cs_beta, X = cb_data$data[[X_dict]]$X, Y = cb_data$data[[j]]$Y, 
                                          XtX = cb_data$data[[X_dict]]$XtX, XtY = cb_data$data[[j]]$XtY, 
                                          YtY = cb_data$data[[j]]$YtY, N = cb_data$data[[j]]$N, 
                                          miss_idx = cb_data$data[[j]]$variable_miss)
                last_profile <- extract_last(cb_obj$cb_model[[j]]$profile_loglike_each)
                change <- abs(cs_profile - last_profile) 
                # - add hoc
                if (min(cb_obj$cb_model[[j]]$multi_correction_univariate[cs_variants]) >= cut){ change <- 0 }
                # ------
                # - check_null
                check_cs_change[i,j] <- change  / diff(range(cb_obj$cb_model[[j]]$profile_loglike_each))
                cs_change[i,j] <- change # / diff(range(cb_obj$cb_model[[j]]$profile_loglike_each))
                first_profile <- extract_first(cb_obj$cb_model[[j]]$profile_loglike_each)
                max_change[i,j] <- (first_profile - last_profile) >= cb_obj$cb_model[[j]]$check_null_max
            } else if (check_null_method == "obj"){
                res <- update_res(X = cb_data$data[[X_dict]]$X, Y = cb_data$data[[j]]$Y, 
                                  XtX = cb_data$data[[X_dict]]$XtX, XtY = cb_data$data[[j]]$XtY, 
                                  N = cb_data$data[[j]]$N, cs_beta,
                                  miss_idx = cb_data$data[[j]]$variable_miss)
                cs_obj <- get_cs_obj(cs_beta, res, cb_obj$cb_model_para$tau, cb_obj$cb_model_para$func_prior, 
                                     cb_obj$cb_model_para$lambda, adj_dep = adj_dep,
                                     LD_obj = cb_obj$cb_model_para$LD_obj,
                                     X = cb_data$data[[X_dict]]$X, N = cb_data$data[[j]]$N, 
                                     XtX = cb_data$data[[X_dict]]$XtX, XtY = cb_data$data[[X_dict]]$XtY, 
                                     YtY = cb_data$data[[X_dict]]$YtY,
                                     miss_idx = cb_data$data[[j]]$variable_miss)
                last_obj <- min(cb_obj$cb_model[[j]]$obj_path)
                change <- abs(cs_obj - last_obj) 
                if (length(cb_obj$cb_model[[j]]$obj_path)==1){
                  total_obj <- 1
                } else {
                  total_obj <- diff(range(cb_obj$cb_model[[j]]$obj_path))
                }
                check_cs_change[i,j] <- change / total_obj
                cs_change[i,j] <- change
                max_change[i,j] <- total_obj >= cb_obj$cb_model[[j]]$check_null_max
            }
        }
    }
    if (!weaker_ucos){ 
        check_cs_change <- cs_change 
        check_null <- cb_obj$cb_model[[j]]$check_null_max
    } 
    cs_change <- as.data.frame(cs_change)
    is_non_null <- which(rowSums( (check_cs_change >= check_null) * max_change ) != 0)
    
    ll = list("cs_change" = cs_change,
              "is_non_null" = is_non_null)
    return(ll)
}


get_purity = function (pos, X=NULL, Xcorr=NULL, N = NULL, n = 100) {
    get_upper_tri = Rfast::upper_tri
    get_median    = Rfast::med
    
    if (sum(is.na(pos))!=0){
        pos <- as.numeric(na.omit(pos))
    }
    
    if (length(pos) == 1)
        return(c(1,1,1))
    else {
        
        # Subsample the columns if necessary.
        if (length(pos) > n)
            pos = sample(pos,n)
        
        if (is.null(Xcorr)) {
            X_sub = X[,pos]
            X_sub = as.matrix(X_sub)
            corr <- suppressWarnings({get_cormat(X_sub)})
            corr[which(is.na(corr))] = 0
            value = abs(get_upper_tri(corr))
        } else {
            Xcorr <- Xcorr # if (!is.null(N)) Xcorr/(N-1) else Xcorr
            value = abs(get_upper_tri(Xcorr[pos,pos]))
        }
        return(c(min(value),
                 sum(value)/length(value),
                 get_median(value)))
    }
}

# ------ Calculate modularity ----------
get_modularity <- function(Weight, B){
    if (dim(Weight)[1] == 1){
        Q <- 0
    } else {
        W_pos <- Weight * (Weight > 0)
        W_neg <- Weight * (Weight < 0)
        N <- dim(Weight)[1]
        K_pos <- colSums(W_pos)
        K_neg <- colSums(W_neg)
        m_pos <- sum(K_pos)
        m_neg <- sum(K_neg)
        m <- m_pos + m_neg
        cate <- B %*% t(B)
        if (m_pos == 0 & m_neg == 0){
            Q <- 0
        } else {
            if (m_pos == 0){
                Q_positive <- 0
                Q_negative <- sum((W_neg - K_neg %*% t(K_neg) / m_neg) * cate) / m_neg
            } else if (m_neg == 0){
                Q_positive <- sum((W_pos - K_pos %*% t(K_pos) / m_pos) * cate) / m_pos
                Q_negative <- 0
            } else {
                Q_positive <- sum((W_pos - K_pos %*% t(K_pos) / m_pos) * cate) / m_pos
                Q_negative <- sum((W_neg - K_neg %*% t(K_neg) / m_neg) * cate) / m_neg
            }
        }
        Q <- m_pos / m * Q_positive - m_neg / m * Q_negative
    }
}


get_n_cluster <- function(hc, Sigma, m=ncol(Sigma), between_cluster = 0.8){
    if (min(Sigma) > between_cluster){
        IND = 1
        Q = 1
    } else {
        Q <- c()
        if (ncol(Sigma) < 10){m = ncol(Sigma)}
        for(i in 1:m)
        {
            index=cutree(hc,i)
            B=sapply(1:i, function(t) as.numeric(index==t))
            Q[i] <- get_modularity(Sigma, B)
        }
        
        IND=which(Q==max(Q))
        L=length(IND)
        if (L>1) IND=IND[1]
    }
    return(list("n_cluster" = IND,
                "Qmodularity" = Q))
}

w_purity <- function(weights, X=NULL, Xcorr=NULL, N = NULL, n = 100, coverage = 0.95, 
                     min_abs_corr = 0.5, median_abs_corr = NULL, miss_idx = NULL){
    
    tmp <- apply(weights, 2, w_cs, coverage = coverage)
    tmp_purity <- apply(tmp, 2, function(tt){
        pos <- which(tt == 1)
        # deal with missing snp here
        if (!is.null(Xcorr)){
            pos <- match(pos, setdiff(1:length(tmp), miss_idx))
        }
        get_purity(pos, X=X, Xcorr=Xcorr, N=N, n = n)
    })
    if (is.null(median_abs_corr)) {
        is_pure = which(tmp_purity[1,] >= min_abs_corr)
    } else {
        is_pure = which(tmp_purity[1,] >= min_abs_corr | tmp_purity[3,] >= median_abs_corr)
    }
    return(is_pure)
}


# - Calculate purity between two confidence sets
get_between_purity = function (pos1, pos2, X=NULL, Xcorr=NULL, N = NULL, miss_idx = NULL, P = NULL){
    
    get_matrix_mult <- function(X_sub1, X_sub2){
        
        X_sub1 = t(X_sub1)
        X_sub2 = t(X_sub2)
        # Standardize each variable
        X_sub1 = X_sub1 / sqrt(rowSums(X_sub1^2))
        X_sub1[is.nan(X_sub1)] = 0
        X_sub2 = X_sub2 / sqrt(rowSums(X_sub2^2))
        X_sub2[is.nan(X_sub2)] = 0
        # Calculate correlations
        cr = tcrossprod(X_sub1, X_sub2)
        return(cr)
        
    }
    get_median = Rfast::med
    
    if (is.null(Xcorr)){
        X_sub1 = X[,pos1,drop=FALSE]
        X_sub2 = X[,pos2,drop=FALSE]
        value <- abs(get_matrix_mult(X_sub1, X_sub2))
        
    } else {
        pos1 <- na.omit(match(pos1, setdiff(1:P, miss_idx)))
        pos2 <- na.omit(match(pos2, setdiff(1:P, miss_idx)))
        # scaling_factor <- if (!is.null(N)) N-1 else 1
        if (length(pos1)!=0 & length(pos2)!=0){
            value = abs(Xcorr[pos1,pos2]) # / scaling_factor
        } else {
            value = 0
        }
        
    }
    return(c(min(value), max(value), get_median(value)))
}






