## BJ statistics
getBJKey <- function(n,alpha,oneSide){
    paste(n,alpha,oneSide)
}

interpolation <- function(x, x1,y1,x2,y2){
    y1+(y2-y1)*(x-x1)/(x2-x1)
}

genericCritical<-function(pvalueFunc,searchRange, n,alpha, oneSide = TRUE){
    rootFunc=function(stat) 
        vapply(stat, function(stat)
            pvalueFunc(stat=stat,n=n, oneSide = oneSide)-alpha,numeric(1))
    res=uniroot(rootFunc,searchRange,extendInt = "yes")
    res$root
}

BJExactCritical<-function(n,alpha, oneSide = TRUE){
    n_l <- floor(n)
    n_h <- n_l+1
    key1 <- paste(n_l,alpha,oneSide)
    key2 <- paste(n_h,alpha,oneSide)
    if(exists(key1,envir = BJCriticalSpace) && exists(key2,envir = BJCriticalSpace)){
        y1 <- get(key1,envir = BJCriticalSpace)
        y2 <- get(key2,envir = BJCriticalSpace)
        res <- interpolation(n,n_l,y1,n_h,y2)
    }else{
        key <- paste(n,alpha,oneSide)
        res <- genericCritical(
        pvalueFunc = BJExactPvalue,searchRange=c(0,1),
        n=n,alpha=alpha,
        oneSide = oneSide
        )
        assign(key,res,envir = BJCriticalSpace)
    }
    res
}

BJExactPvalue<-function(statValue,n, oneSide = TRUE){
    bounds <- BJLocalCritical(statValue=statValue,n=n)
    if(oneSide){
        bounds$h <- rep(1,length(bounds$h))
    }
    1-orderedProb(bounds$l, bounds$h)
}


updataBJStat <- function(statValue, samples, x){
    x[x==0] <- 1e-10
    x[x==1] <- 1-1e-10
    samples <- sort(c(samples,x))
    i <- seq_along(samples)
    -min(pbeta(samples, i, length(samples)-i+1))
}

BJCritical <- function(alpha, n){
    if(n<=1000){
        -BJExactCritical(n,alpha)
    }else{
        log(1-alpha)/(2*log(n)*log(log(n)))
    }
}

BJGlobal <- function(){
    .GlobalNullStat(statUpdataFunc = updataBJStat,
                    criticalFunc = BJCritical)
}

