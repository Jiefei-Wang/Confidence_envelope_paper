test_that(
    "fast algorithm",
    {
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

        statFunc <- function(pvalues){
            sum(qnorm(pvalues))
        }

        statCriticalFunc <- function(n, alpha){
            qnorm(alpha)*sqrt(n)
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

        I <- rep(list(NULL), groups)
        for(i in seq_along(I)){
            I[[i]] <- (groupSize*(i-1)+1):(groupSize*i)
        }

        for(i in 1:10){
            pvalues <- runif(n)
            pvalues[trueAlter] <- rbeta(length(trueAlter), 0.5, 2)

            nullSets <- findNullSets(pvalues, mytest, alpha)

            k<- c()
            for(i in 1:groups){
                k[i] <- maxLoss(lossFuncSets[[i]], nullSets)
            }

            param <- GCS_FP_param(pvalues, I,statFunc)
            k_candidate<- expand.grid(1:groupSize, 1:groupSize,1:groupSize,1:groupSize)

            for(i in 1:nrow(k_candidate)){
                lossFunc <- lossGenerator2(k_candidate[i,], groupCounts, groupSize)
                loss1 <- maxLoss(lossFunc, nullSets)
                loss2 <- GCS_FP(param, k_candidate[i,], statFunc,statCriticalFunc, alpha)
                expect_equal(loss1,loss2)
            }
        }
    }
)
