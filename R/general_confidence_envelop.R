#' Compute the confidence envelop for any function of the true null hypotheses
#' 
#' @param pvals A named numeric vector of p-values for the testing of the null hypotheses
#' @param func A function that takes a vector of indices of the true null hypotheses and returns a numeric value
#' @param alpha The significance level
#' @param nullTest A function that takes a vector of p-values and alpha level and returns significance of the global test
generalCE <- function(pvals, func, alpha = 0.05, nullTest = stoufferTest){
    hypos <- names(pvals)
    if (is.null(hypos)) hypos <- seq_along(pvals)
    lower <- func(integer(0))
    upper <- func(integer(0))
    nonSigSetLower <- integer(0)
    nonSigSetUpper <- integer(0)
    n <- length(hypos)
    for(k in rev(seq_along(hypos))){
        for (i in seq_len(choose(n, k))){
            if (i == 1)
                cbn <- 1:k
            else
                cbn <- gen.next.cbn(cbn, n)
            hyposSubset <- hypos[cbn]
            psubset <- pvals[hyposSubset]
            globalPval <- nullTest(psubset, alpha)
            if (!globalPval) {
                value <- func(cbn)
                if (value < lower) {
                    lower <- value
                    nonSigSetLower <- cbn
                }
                if (value > upper) {
                    upper <- value
                    nonSigSetUpper <- cbn
                }
            }
        }
    }
    list(
        lower = lower,
        upper = upper,
        nonSigSetLower = nonSigSetLower,
        nonSigSetUpper = nonSigSetUpper
    )
}

generalCE2 <- function(pvals, func, alpha = 0.05, nullTest = stoufferTest){
    hypos <- names(pvals)
    if (is.null(hypos)) hypos <- seq_along(pvals)
    if(all(pvals<= 1-alpha/length(pvals))){
        lower <- func(integer(0))
        upper <- func(integer(0))
        nonSigSetLower <- integer(0)
        nonSigSetUpper <- integer(0)
    }else{
        message("hit")
        lower <- Inf
        upper <- -Inf
        nonSigSetLower <- NA
        nonSigSetUpper <- NA
    }
    n <- length(hypos)
    for(k in rev(seq_along(hypos))){
        for (i in seq_len(choose(n, k))){
            if (i == 1)
                cbn <- 1:k
            else
                cbn <- gen.next.cbn(cbn, n)
            hyposSubset <- hypos[cbn]
            psubset <- pvals[hyposSubset]
            globalPval <- nullTest(psubset, alpha)
            if (!globalPval) {
                value <- func(cbn)
                if (value < lower) {
                    lower <- value
                    nonSigSetLower <- cbn
                }
                if (value > upper) {
                    upper <- value
                    nonSigSetUpper <- cbn
                }
            }
        }
    }
    list(
        lower = lower,
        upper = upper,
        nonSigSetLower = nonSigSetLower,
        nonSigSetUpper = nonSigSetUpper
    )
}


#' stouffer test    
stoufferTest <- function(pvalues, alpha){
    poolr::stouffer(runif(10))$p
}

#' Bonferroni test
bonferroniTest <- function(pvalues, alpha){
    sum(pvalues < alpha/length(pvalues)) > 0
}

BJTestGenerator <- function(l=c(0,1), h=NULL){
    function(pvals, alpha){
        critical <- BJCritical(alpha, length(pvals))
        statValue <- BJStat(pvals, l, h)
        statValue < critical
    }
}


pi_fast <- function(pvals, alpha, lambda){

}

#' @export 
pi_BJ_fast <- function(func, pvals, alpha, l=c(0,1), h=NULL){
    lower <- func(integer(0))
    upper <- func(integer(0))
    nonSigSetLower <- integer(0)
    nonSigSetUpper <- integer(0)
    pvals_sorted <- sort(pvals, decreasing = TRUE)
    ## the original index of the sorted p-values
    idx <- order(pvals, decreasing = TRUE)
    ## hypotheses of set k
    for(k in seq_along(pvals)){
        critical <- BJCritical(alpha, k, l=l,h=h)
        ## The bound of the p-values
        indexL <- getIndex(k,l)
        indexH <- getIndex(k,h)
        bounds <- BJLocalBounds(critical, k)
        lb <- bounds$l
        hb <- bounds$h
        lb[setdiff(seq_len(k),indexL)] <- 0
        hb[setdiff(seq_len(k),indexH)] <- 1
        lb <- rev(lb)
        hb <- rev(hb)

        pselected <- rep(NA, k)
        select_i <- 1
        for(candidate_i in seq_along(pvals_sorted)){
            cur_l <- lb[select_i]
            cur_h <- hb[select_i]
            cur_p <- pvals_sorted[candidate_i]
            if(cur_p>=cur_l && cur_p<=cur_h){
                pselected[select_i] <- candidate_i
                select_i = select_i + 1
            }
            if(select_i > k){
                break
            }
        }
        ## cannot find a p-value in the range, go to next k
        if (any(is.na(pselected)))
            next

        cbn <- idx[pselected]
        value <- func(cbn)
        if (value > upper) {
            upper <- value
            nonSigSetUpper <- cbn
        }
    }
    list(
        lower = lower,
        upper = upper,
        nonSigSetLower = nonSigSetLower,
        nonSigSetUpper = nonSigSetUpper
    )
}
