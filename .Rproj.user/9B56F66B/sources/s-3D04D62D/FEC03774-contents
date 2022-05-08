library(doParallel)
library(foreach)
source("R/exceedance.R")
registerDoParallel(cores = 24)

alpha <- 0.2

n <- 12
trueAlter <- 1:12
trueNull <- setdiff(1:n, trueAlter)
groups <- 4
groupSize <- n/groups

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
mytest <- function(pvalues){
    library("poolr")
    stouffer(pvalues)$p
}


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

simResult <- c()
simResult <- foreach(i = 1:100, .combine = rbind) %dopar%{
    # for(i in 1:1000){

    pvalues <- runif(n)
    pvalues[trueAlter] <- rbeta(length(trueAlter), 0.5, 2)

    nullSets <- findNullSets(pvalues, mytest, alpha)

    result<- c()
    for(i in 1:groups){
        result[i] <- maxLoss(lossFuncSets[[i]], nullSets)
    }

    k <- result + 1L
    tmp_k <- k
    loss <- 0L
    while(loss<=1){
        for(i in 1:groups){
            tmp_k <- k
            tmp_k[i] <- tmp_k[i] - 1
            lossFunc_sharper <- lossGenerator2(tmp_k, groupCounts, groupSize)
            loss <- maxLoss(lossFunc_sharper, nullSets)
            if(loss <= 1){
                k <- tmp_k
            }
        }
    }

    c(result, k - 1)
}

colMeans(simResult)


generalLoss <- simResult[,1]
individualLoss <- simResult[,2:ncol(simResult)]
## Probability of the bound
mean(trueLoss1 <= simResult[,1])
mean(trueLoss2 <= simResult[,2])
mean(trueLoss12 <= simResult[,3])



which(generalLoss < trueLoss & simResult[,3]==1)

table(generalLoss, simResult[,3])





mean(generalLoss)


mean(rowSums(individualLoss) < trueLoss)

sapply(nullSets,lossFuncAll)




