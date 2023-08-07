library(foreach)
library(doParallel)
library(parallel)
library(doRNG)
registerDoParallel(cores = detectCores())

## Use BJ statistic as the local test method
globalTest <- BJGlobal()

## Significance level
alpha <- 0.05

## Number of hypotheses
n <- 500
## Number of groups
nGroups <- 5
groupSize <- n/nGroups
## Number of true alternatives in each group
nGroupAlter <- 50
nGroupNull <- groupSize - nGroupAlter

## Index of hypothesis in each group
groupsIndex <- lapply(1:nGroups, function(i) 1:groupSize + (i-1)*groupSize)

## Indice of the true alternative hypotheses
trueAltIndex <- rep(nGroupNull + seq_len(nGroupAlter), nGroups) + rep((0:(nGroups-1))*groupSize,each=nGroupAlter)

## Number of simulations
nSim <- 2
## Alternative distribution(Beta distribution)
a <- 0.5
b <- 1

## False coverage number you want to control
FC <- 1
registerDoRNG(seed = 123)
final <- foreach(i = 1:nSim, .combine = rbind) %dopar% {
    devtools::load_all()
    pvalues <- runif(n)
    pvalues[trueAltIndex] <- rbeta(length(trueAltIndex), a,b)
    res <- findOptimalBounds(globalTest, pvalues, groupsIndex, alpha, FC = FC)
    res$cutoff
}

## Average of the sum of the true discovery guarantee
mean(rowSums(final))

