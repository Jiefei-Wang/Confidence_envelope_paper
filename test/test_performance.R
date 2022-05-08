library(doParallel)
library(foreach)
source("R/exceedance.R")
registerDoParallel(cores = 24)

alpha <- 0.2

n <- 16
trueAlter <- 1:16
trueNull <- setdiff(1:n, trueAlter)
groups <- 4
groupSize <- n/groups



statFunc <- function(pvalues){
    sum(qnorm(pvalues))
}

statCriticalFunc <- function(n, alpha){
    qnorm(alpha)*sqrt(n)
}


groupCountGenerator <- function(){
    groups <- groups
    groupSize <- groupSize
    function(nullHypo){
        counts <- rep(0L, groups)
        ind <- (nullHypo-1L)%/%groupSize + 1L
        for(i in ind){
            counts[i] <- counts[i] + 1L
        }
        counts
    }
}

groupCounts <- groupCountGenerator()





lossGenerator1 <- function(i){
    force(i)
    groupCounts <- groupCounts
    groupSize <- groupSize
    # fdpCritical <- fdpCritical
    function(nullHypo){
        counts <- groupCounts(nullHypo)
        sum(counts[i])
    }
}

lossGenerator2 <- function(k, groupCounts, groupSize){
    force(k)
    function(nullHypo){
        counts <- groupCounts(nullHypo)
        sum(counts >= k)
    }
}
#
lossFuncSets <- list()
for(i in 1:groups){
    lossFuncSets[[i]] <- lossGenerator1(i)
}


trueLoss1 <- lossFuncSets[[1]](trueNull)
trueLoss2 <- lossFuncSets[[2]](trueNull)

I <- rep(list(NULL), groups)
for(i in seq_along(I)){
    I[[i]] <- (groupSize*(i-1)+1):(groupSize*i)
}


simResult <- c()
simResult <- foreach(i = 1:100, .combine = rbind) %dopar%{
    # for(i in 1:1000){

    pvalues <- runif(n)
    pvalues[trueAlter] <- rbeta(length(trueAlter), 0.5, 2)



    param <- GCS_FP_param(pvalues, I, statFunc)
    k_candidate<- expand.grid(1:groupSize, 1:groupSize,1:groupSize,1:groupSize)

    system.time(
        {
            nullSets <- findNullSets(pvalues, mytest, alpha)
            for(i in 1:nrow(k_candidate)){
                lossFunc <- lossGenerator2(k_candidate[i,], groupCounts, groupSize)
                loss1 <- maxLoss(lossFunc, nullSets)
            }
        }
    )

    system.time(
        {
            for(i in 1:nrow(k_candidate)){
                loss2 <- GCS_FP(param, k_candidate[i,], statFunc, statCriticalFunc, alpha)
            }

        }
    )

    c(result, k - 1)
}

colMeans(simResult)
