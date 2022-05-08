
## Reject: pvalue <= alpha
## not reject: pvalue > alpha

#' Computing the confidence envolop of the false discover proportion for the data
#' 
#' Computing the 1 - alpha level confidence envolop of the false discover 
#' proportion(FDP) given a set of rejected hypotheses. 
#' The confidence envolop can be viewed as a measurement of the quality of 
#' the statistical inference.
#' 
#' @param profiled_data an exceedance_profile object
#' @param alpha numeric, the confidence level
#' @param ri integer, the index of the rejected hypotheses, see details.
#' @param sri integer, the index of the ascending ordered p-values which the 
#' corresponding hypotheses are rejected, see details.
#' @param rx numeric, the value of the pvalues which the 
#' corresponding hypotheses are rejected, see details.
#' 
#' @details 
#' This function is for constructing the confidence envolop of the
#' FDP given the set of rejected hypothese. The confidence envolop
#' depends on three factors: 
#' \itemize{
#' \item The p-value samples
#' \item The confidence level `alpha`
#' \item The rejected hypotheses.
#' }
#' Therefore, given the data, confidence level and the hypotheses that you want
#' to reject, we can obtain a `1 - alpha` confidence envolop of the FDP. 
#' 
#' The rejected hypotheses can be expressed in three ways. You can use the 
#' original index `ri` to indicate which hypotheses you want to reject. For
#' example, if `ri = 1:2`, it means the first and second hypotheses are rejected.
#' 
#' However, in practice, it is more common to reject the hypotheses which 
#' have small pvalues. You can achieve it by providing the parameter `sri`. 
#' For example, if `sri = 1:2`, it means the hypothese which have the smallest 
#' or second smallest pvalues are rejected. Alternatively, `rx` can be used if
#' you want to match the pvalues not the index. That is, a hypotheis is 
#' rejected if its pvalue matches any value in `rx`.
#' 
#' @return a `1 - alpha` level confidence envolop
#' @inherit exceedance_inference examples
#' @export
exceedance_confidence<-function(profiled_data, alpha,ri=NULL,sri = NULL,rx=NULL){
    x_rank <- profiled_data$profile$x_rank
    sorted_i <- as.integer(get_ordered_index(x_rank,ri,sri,rx))
    
    
    confidence_func <- get(profiled_data$param$confidence_func)
    
    result <- confidence_func(
        profiled_data = profiled_data, alpha = alpha, 
        sorted_i = sorted_i)
    result
    
}

confidence_general_GW_general <- function(profiled_data,alpha,
                                          sorted_i){
    profile <- profiled_data$profile
    param <- profiled_data$param
    m <- profile$m
    unreject_set<- profile$unreject_set
    
    FP <- general_GW_compute_FP(unreject_set,m,sorted_i, alpha)
    FDR <- FP/length(sorted_i)
    FDR
}

confidence_general_GW_JW <- function(profiled_data,alpha,
                                     sorted_i){
    profile <- profiled_data$profile
    param <- profiled_data$param
    pvalue_func <- param$pvalue_func
    m <- profile$m
    x_sort<- profile$x_sort
    cache <- profile$cache
    
    candidate_set <- setdiff(seq_len(m),sorted_i)
    max_candidate_num <- length(candidate_set)
    total_positive <- length(sorted_i)
    
    # sorted_i_key <- digest(sorted_i)
    
    FDR <- 0
    ## pick false positive number
    for(FP in rev(seq_len(total_positive))){
        ## pick the size of the candidate set(exclude the false positive)
        for(unreject_size in seq_len(max_candidate_num+1L)-1L){
            ## construct the entire candidate set
            if(unreject_size==0L){
                index <- sorted_i[seq(length(sorted_i)-FP + 1L,length(sorted_i))]
            }else{
                index <- c(sorted_i[seq(length(sorted_i)-FP + 1L,length(sorted_i))],
                           candidate_set[seq(max_candidate_num-unreject_size + 1L,max_candidate_num)]
                )
            }
            index <- sort(index)
            pvalue <- pvalue_func(x_sort[index])
            ## test if the candidate set is rejected at level alpha
            if(pvalue > alpha){
                FDR <- FP/length(sorted_i)
                break;
            }
        }
        if(pvalue > alpha){
            break;
        }
    }
    FDR
}


## Need optimization
confidence_fast_GW_kth_p_index <- function(profiled_data,alpha,
                                           sorted_i){
    profile <- profiled_data$profile
    param <- profiled_data$param
    max_alpha <- profile$max_alpha
    k<- param$param1
    m <- profile$m
    local_level <- profile$local_level
    
    
    if(alpha>max_alpha)
        return(min(k-1,length(sorted_i))/length(sorted_i))
    
    J_all <- which(local_level[k:m]>alpha)
    
    if(length(J_all)==0){
        J <- m + 1L
    }else{
        J <- min(J_all)+k-1
    }
    
    FP <- sum(sorted_i<k)+sum(sorted_i>=J)
    # U <- setdiff(seq_len(m),k-1+seq_len(J-k))
    FDR <- FP/length(sorted_i)
    FDR
}

confidence_fast_GW_kth_p_proportion <- function(profiled_data,alpha,
                                                sorted_i){
    profile <- profiled_data$profile
    param <- profiled_data$param
    
    q<- param$param1
    m <- profile$m
    x_sort <- profile$x_sort
    
    k_list <- pmax(ceiling(seq_len(m)*q),1L)
    cut <- qbeta(alpha,k_list,seq_len(m)-k_list+1L)
    # browser()
    nonsig_index <- rep(m+1L,m)
    for(i in seq_along(nonsig_index)){
        ind <- which(x_sort>cut[i])
        if(length(ind!=0)){
            if(ind[1]<=m-i+k_list[i]){
                nonsig_index[i] <- max(ind[1],k_list[i])
            }
        }
    }
    # browser()
    FP <- 0
    for(n in seq_along(cut)){
        cur_cut <- cut[n]
        cur_k <- k_list[n]
        index <- nonsig_index[n]
        ## all sets are rejected
        if(index == m+1L)
            next
        L <- sum(sorted_i<index)
        if(L<=cur_k-1){
            FP <- length(sorted_i)
            break
        }else{
            if(L!=n){
                FP <- max(FP,cur_k-1L+length(sorted_i)-L)
            }else{
                FP <- max(FP,cur_k-1)
            }
        }
    }
    FP/length(sorted_i)
}







confidence_fast_GW_order_general <- function(profiled_data,alpha,
                                             sorted_i){
    profile <- profiled_data$profile
    param <- profiled_data$param
    range_type <- param$range_type
    param1 <- param$param1
    param2 <- param$param2
    statistic <- param$statistic
    m <- profile$m
    x_sort<- profile$x_sort
    params_key <- profile$params_key
    
    preprocessed_key <- paste0("GW",alpha,params_key)
    # browser()
    rj_num <- length(sorted_i)
    ## reduce the loop number by checking the current FDR
    FDR <- 0
    for(n in rev(seq_len(m))){
        # message(n)
        if(range_type=="proportion"){
            indexL <- get_index_from_proportion(n=n,param=param1)
            indexU <- get_index_from_proportion(n=n,param=param2)
        }else{
            indexL <- param1[param1<=n]
            indexU <- param2[param2<=n]
            ## Check if the current sample size satisfies
            ## The minimum requirement of the test.
            if(length(indexL)==0&&length(indexU)==0){
                FDR <- max(FDR, min(n,rj_num) / rj_num)
                next
            }
        }
        
        critical <- get_critical(statName= statistic, n=n, alpha=alpha, 
                                 indexL=indexL, indexU=indexU)
        bound <- get_local_critical(statName = statistic, n= n, critical=critical,
                                    indexL=indexL,indexU=indexU)
        
        # x_range <- get_range_by_bound(sx=x_sort,bound=bound)
        x_range <- C_get_range_by_bound(R_sx=x_sort,R_l=bound$l,R_h=bound$h)
        if(is.null(x_range)){
            next
        }
        P <- x_range$P
        Q <- x_range$Q
        # browser()
        FDR <- max(FDR,C_GW_compute_FDR(sorted_i,P,Q,rj_num,n))
        
    }
    FDR
}

confidence_combine_GW <- function(profiled_data,alpha,
                                  sorted_i){
    profiles <- profiled_data$profile$profiles
    alpha_weight <- profiled_data$param$alpha_weight
    
    alphas <- alpha/sum(alpha_weight)*alpha_weight
    confidences <- lapply(seq_along(profiles), 
                          function(i,profiles,alphas)
                              exceedance_confidence(profiled_data=profiles[[i]],
                                                    alpha=alphas[i],
                                                    sri = sorted_i),
                          profiles=profiles,
                          alphas=alphas)
    min(as.numeric(confidences))
}