## BJ statistics
getBJKey <- function(n,alpha,l,h){
  paste(n, alpha, paste0(l,collapse = ","), ";",paste0(h,collapse = ","))
}
getBJLocalKey <- function(statValue, n){
  paste0(statValue,";",n)
}

interpolation <- function(x, x1,y1,x2,y2){
    y1+(y2-y1)*(x-x1)/(x2-x1)
}

genericCritical<-function(pvalueFunc,searchRange, n,alpha, l, h){
    rootFunc=function(stat) 
        vapply(stat, function(stat)
            pvalueFunc(stat=stat,n=n, l = l, h = h)-alpha,numeric(1))
    res=uniroot(rootFunc,searchRange,extendInt = "yes")
    res$root
}

BJPlusLocal <- function(pvals_sorted) {
    n <- length(pvals_sorted)
    vapply(seq_len(n), function(i)
        pbeta(pvals_sorted[i], i, n - i + 1), numeric(1))
}
BJMinusLocal <- function(pvals_sorted) {
    1 - BJPlusLocal(pvals_sorted)
}
BJStat <- function(pvals, l = c(0,1), h=NULL) {
    indexL <- getIndex(n,l)
    indexH <- getIndex(n,h)
    pvals_sorted <- sort(pvals)
    BJPlus <- BJPlusLocal(pvals_sorted)
    BJMinus <- BJMinusLocal(pvals_sorted)
    min(c(BJPlus[indexL],BJMinus[indexH]), na.rm = TRUE)
}



BJExactCritical<-function(n,alpha, l = c(0,1), h = NULL){
    key <- getBJKey(n, alpha,l,h)
      
    if(exists(key,envir = BJCriticalSpace)){
        res <- get(key,envir = BJCriticalSpace)
    }else{
        res <-genericCritical(
            pvalueFunc = BJExactPvalue,
            searchRange=c(0,1),
            n=n,
            alpha=alpha,
            l = l,
            h = h
            )
        assign(key,res,envir = BJCriticalSpace)
    }
    res
}

getIndex <- function(n, l){
    if(length(l)==0){
        return(integer(0))
    }
    stopifnot(length(l)==2)
    l_range <- l*(n-1)+1
    idx <- seq_len(n)
    idx[idx >= l_range[1] & idx <= l_range[2]]
}

BJExactPvalue <- function(statValue, n, l = c(0,1), h = NULL){
    indexL <- getIndex(n,l)
    indexH <- getIndex(n,h)
    bounds <- BJLocalBounds(statValue=statValue,n=n)
    bounds$l[setdiff(seq_len(n),indexL)] <- 0
    bounds$h[setdiff(seq_len(n),indexH)] <- 1
    1-orderedProb(bounds$l, bounds$h)
}

## The boundary of the ordered pvalues
BJLocalBounds<-function(statValue, n){
    key <- getBJLocalKey(statValue, n)
    
    if(exists(key,envir = BJCriticalSpace)){
      res <- get(key,envir = BJCriticalSpace)
    }else{
      l=vapply(seq_len(n),function(x)qbeta(statValue,x,n-x+1),numeric(1))
      h=vapply(seq_len(n),function(x)qbeta(1 - statValue,x,n-x+1),numeric(1))
      res <- list(l =l,h= h)
      assign(key,res,envir = BJCriticalSpace)
    }
    res
}

## The boundary of the ordered pvalues at the critical value
## one-sided only
BJLocalCriticals <- function(n, alpha, l = c(0,1), h = NULL){
    critical <- BJCritical(alpha, n, l=l, h=h)
    
    indexL <- getIndex(n,l)
    indexH <- getIndex(n,h)
    bounds <- BJLocalBounds(critical, n)
    bounds$l[setdiff(seq_len(n),indexL)] <- 0
    bounds$h[setdiff(seq_len(n),indexH)] <- 1
    bounds
}


BJCritical <- function(alpha, n, l = c(0,1), h = NULL){
    if(n>1000&&identical(l,c(0,1))&&is.null(h)){
        log(1-alpha)/(2*log(n)*log(log(n)))
    }else{
        BJExactCritical(n,alpha, l = l, h = h)
    }
}

BJGlobal <- function(){
    .GlobalNullStat(statUpdataFunc = updataBJStat,
                    criticalFunc = BJCritical)
}

