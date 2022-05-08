#' Compute the critical value of generalized Kolmogorov-Smirnov test statistics
#' 
#' Compute the critical value of generalized Kolmogorov-Smirnov test statistics,
#' see details on how to use the function.
#' 
#' @param alpha numeric, the type I error rate for the critical value. Please do 
#' not be confused with `alpha0`.
#' @param n Integer, the sample size of the data.
#' 
#' @examples 
#' ## Compute the critical value of the KS test
#' ## of sample size 10
#' GKSCritical(alpha = 0.05, n = 10, statName = "KS")
#' 
#' ## The critical value for the test that
#' ## only considers the first 3 ordered samples
#' ## All gives the same result.
#' GKSCritical(alpha = 0.05, n = 10, alpha0 = 0.3, statName = "KS")
#' GKSCritical(alpha = 0.05, n = 10, index = 1:3, statName = "KS")
#' GKSCritical(alpha = 0.05, n = 10, indexL = 1:3, indexU = 1:3, statName = "KS")
#' 
#' 
#' @return A critical value
#' @inheritParams GKSStat
#' @inherit GKSStat details
#' @rdname critical
#' @export
GKSCritical <-function(n,alpha,
                       index= NULL, 
                       indexL=NULL,indexU= NULL,
                       statName = c("KS","KS+","KS-","BJ","BJ+","BJ-","HC","HC+","HC-","Simes")){
    statName <- match.arg(statName)
    stopifnot(!is.null(alpha))
    stopifnot(!is.null(n))
    
    idx <- getGKSIndex(statName = statName, n = n, 
                       index= index, 
                       indexL = indexL, indexU = indexU)
    indexL <- idx$indexL
    indexU <- idx$indexU
    statName <- idx$statName
    
    statValue <- call_func(root = "Critical",prefix = statName,
                           n=n,
                           alpha=alpha, 
                           indexL=indexL,
                           indexU=indexU)
    return(statValue)
}


genericCritical<-function(pvalueFunc, searchRange,
                          n,alpha,
                          indexL,indexU){
    rootFunc=function(stat) 
        vapply(stat, function(stat)
            pvalueFunc(stat=stat,n=n,indexL=indexL,indexU=indexU)-alpha,numeric(1))
    res=uniroot(rootFunc,searchRange,extendInt = "yes")
    res$root
}


HCCritical<-function(n,alpha,indexL=seq_len(n),indexU=seq_len(n)){
    genericCritical(
        pvalueFunc= HCPvalue,searchRange=c(0,100),
        n=n,alpha=alpha,
        indexL=indexL,indexU= indexU
    )
}

BJCritical<-function(n,alpha,indexL=seq_len(n),indexU=seq_len(n)){
    genericCritical(
        pvalueFunc= BJPvalue,searchRange=c(0,1),
        n=n,alpha=alpha,
        indexL=indexL,indexU= indexU
    )
}

KSCritical<-function(n,alpha,indexL=seq_len(n),indexU=seq_len(n)){
    genericCritical(
        pvalueFunc= KSPvalue,searchRange=c(0,1),
        n=n,alpha=alpha,
        indexL=indexL,indexU= indexU
    )
}

SimesCritical<-function(n,alpha,indexL=seq_len(n),indexU=seq_len(n)){
    alpha
}