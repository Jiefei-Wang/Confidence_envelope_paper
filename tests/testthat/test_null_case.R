nsims <- 1000
n <- 6
n_null <- 6
a <- 1
b <- 10
trueNulls <- seq_len(n_null)

generatePvals <- function(n, n_null, a, b){
    pval_null <- runif(n_null)
    pval_alt <- rbeta(n - n_null, a,b)
    c(pval_null, pval_alt)
}

generateDepPvals <- function(n, n_null, a, b){
    pval_null <- rep(runif(1), n_null)
    pval_alt <- rep(rbeta(1, a,b), n - n_null)
    c(pval_null, pval_alt)
}

positive_proportion <- function(nHypos, selection){
    function(nulls){
        ## if no true alternative hypothesis, return 1
        if (length(nulls) == nHypos) return(1)

        (length(selection)-sum(nulls%in%selection))/(nHypos - length(nulls))
    }
}

funcList <- list(
    positive_proportion(n, 1:3),
    positive_proportion(n, 1:4),
    positive_proportion(n, 4:6),
    positive_proportion(n, 3:6)
)

test_that(
    "All nulls",
    {
        for(func in funcList){
            result <- list()
            for(i in 1:nsims){
                pvals <- generatePvals(n, n_null, a, b)
                result[[i]] <- generalCE2(pvals, func = func, alpha = 0.05, nullTest = bonferroniTest)
            }
            result <- do.call(rbind, result)
            true_coverage <- sum(result[,1] <= func(trueNulls) & result[,2] >= func(trueNulls))/nsims

            testthat::expect_true(true_coverage > 0.95)
        }
    }
)


test_that(
    "All nulls, perfect dependence",
    {
        for(func in funcList){
            result <- list()
            for(i in 1:nsims){
                pvals <- generateDepPvals(n, n_null, a, b)
                result[[i]] <- generalCE2(pvals, func = func, alpha = 0.05, nullTest = bonferroniTest)
            }
            result <- do.call(rbind, result)
            true_coverage <- sum(result[,1] <= func(trueNulls) & result[,2] >= func(trueNulls))/nsims

            testthat::expect_true(true_coverage > 0.95)
        }
    }
)
