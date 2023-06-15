n <- 12
nGroups <- 3
groupSize <- n/nGroups
nGroupAlter <- 2
alpha <- 0.05
globalTest <- FisherGlobal()

trueNullList <- rep(TRUE, n)
if (nGroupAlter > 0) {
    for(k in (1:nGroups) * groupSize){
        trueNullList[k - (1:nGroupAlter) +1] <- FALSE
    }
}

trueNullIndex <- which(trueNullList)
trueAltIndex <- which(!trueNullList)
hypoGroup <- rep(1:nGroups, each = groupSize)

## Given a set of hypotheses, 
## return the number of hypotheses in each group
hypothesisInGroups <- function(hypotheses, AllowedGroups=NULL){
    counts <- rep(0, nGroups)
    for(i in hypotheses){
        counts[hypoGroup[i]] <- counts[hypoGroup[i]] + 1
    }
    if(!is.null(AllowedGroups))
        counts <- counts[AllowedGroups]
    counts
}

lossGenerator <- function(cutoff, AllowedGroups=NULL){
    force(cutoff)
    function(trueNulls){
        if(is.null(AllowedGroups))
            AllowedGroups <- seq_along(cutoff)
        ## Number of true nulls in each group
        counts <- hypothesisInGroups(trueNulls, AllowedGroups)
        ## Whether the number of true alternatives in each group 
        ## is less than cutoff
        sum(groupSize - counts < cutoff[AllowedGroups])
    }
}

        

test_that(
    "fast algorithm, no filter",
    {
        cutoff <- c(4, 2, 1)
        loss <- lossGenerator(cutoff)

        set.seed(1)
        nSim <- 20
        final <- c()
        for(i in 1:nSim){
            pvalues <- rbeta(n, 0.5, 4)
            nullSets <- findNullSets(pvalues, globalTest, alpha)
            res1 <- GE(pvalues, alpha, globalTest, loss, nullSets = nullSets)$max
            res2 <- FalseCoverage(globalTest, pvalues, groups, cutoff, alpha)
            final <- rbind(final, c(res1, res2))
            message(i)
        }
        testthat::expect_identical(final[,1], final[,2])
    }
)


test_that(
    "fast algorithm with filter",
    {
        cutoff <- c(1, 0, 0)
        loss <- lossGenerator(cutoff)
        
        set.seed(1)
        nSim <- 20
        final <- c()
        for(i in 1:nSim){
            pvalues <- rbeta(n, 0.5, 4)
            nullSets <- findNullSets(pvalues, globalTest, alpha)
            res1 <- GE(pvalues, alpha, globalTest, loss, nullSets = nullSets)$max
            res2 <- FalseCoverage(globalTest, pvalues, groups, cutoff, alpha)
            final <- rbind(final, c(res1, res2))
            message(i)
        }
        testthat::expect_identical(final[,1], final[,2])
    }
)
