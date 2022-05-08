library(doParallel)
library(foreach)
source("R/exceedance.R")
registerDoParallel(cores = 24)

alpha <- 0.2

# n <- 10
# trueAlter <- c(1,3,5,7,9)
# trueNull <- setdiff(1:n, trueAlter)
# groups <- 5
# groupSize <- n/groups
# fdpCritical <- 0.1

n <- 6
trueAlter <- 5:6
trueNull <- setdiff(1:n, trueAlter)
groups <- 2
groupSize <- n/groups

groupCountGenerator <- function(){
    groups <- groups
    groupSize <- groupSize
    function(nullHypo){
        counts <- rep(0, groups)
        for(i in nullHypo){
            ind <- (i-1)%/%groupSize + 1
            counts[ind] <- counts[ind] + 1
        }
        counts
    }
}

groupCounts <- groupCountGenerator()
mytest <- function(pvalues){
    library("poolr")
    stouffer(pvalues)$p
}
lossGenerator <- function(i){
    force(i)
    groupCounts <- groupCounts
    groupSize <- groupSize
    # fdpCritical <- fdpCritical
    function(nullHypo){
        counts <- groupCounts(nullHypo)
        sum(counts[i]==groupSize)
    }
}

lossFuncSets <- list()
for(i in 1:groups){
    lossFuncSets[[i]] <- lossGenerator(i)
}
lossFuncAll <- lossGenerator(1:groups)

trueLoss1 <- lossFuncSets[[1]](trueNull)
trueLoss2 <- lossFuncSets[[2]](trueNull)
trueLoss12 <- lossFuncAll(trueNull)
simResult <- c()
simResult <- foreach(i = 1:1000, .combine = rbind) %dopar%{
# for(i in 1:1000){

    pvalues <- runif(n)
    pvalues[trueAlter] <- rbeta(length(trueAlter), 0.5, 2)

    nullSets <- findNullSets(pvalues, mytest, alpha)

    result1 <- maxLoss(lossFuncSets[[1]], nullSets)
    result2 <- maxLoss(lossFuncSets[[2]], nullSets)
    result12 <- maxLoss(lossFuncAll, nullSets)

    rj1 <- result1==0 || result12==0
    rj2 <- result2==0 || result12==0

    rjBon1 <- mytest(pvalues[1:3]) <alpha/2
    rjBon2 <- mytest(pvalues[4:6]) <alpha/2

    c(result1, result2, result12, rj1, rj2,rjBon1,rjBon2)
    # simResult <- rbind(simResult, c(result1, result2,turehit, rj))
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




