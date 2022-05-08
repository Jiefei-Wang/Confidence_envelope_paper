library(doParallel)
library(foreach)
source("R/exceedance.R")
registerDoParallel(cores = 24)

alpha <- 0.2

n <- 6
trueAlter <- 5:6
trueNull <- setdiff(1:n, trueAlter)
groups <- 2
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

simResult <- foreach(i = 1:1000, .combine = rbind) %do%{
    # for(i in 1:100){

    pvalues <- runif(n)
    pvalues[trueAlter] <- rbeta(length(trueAlter), 0.5, 2)

    k <- rep(groupSize, 2)
    param <- GCS_FP_param(pvalues, I, statFunc)

    k1 <- c(groupSize, groupSize+1)
    loss1 <- GCS_FP(param, k1, statFunc, statCriticalFunc, alpha)
    k2 <- c(groupSize+1, groupSize)
    loss2 <- GCS_FP(param, k2, statFunc, statCriticalFunc, alpha)
    k12 <- c(groupSize, groupSize)
    loss12 <- GCS_FP(param, k12, statFunc, statCriticalFunc, alpha)


    c(loss1, loss2, loss12)
}



mean(simResult[,3] <=0)
mean(simResult[,3] <=1)
colMeans(simResult)


test <- data.frame(test1= simResult[,1]>=4, test2=simResult[,2]>=4)

table(test)
