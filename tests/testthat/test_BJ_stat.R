
l <- c(0,0.5)
h <- c(0.5,1)

test_that(
    "BJ statistics and pvalue, l=c(0,1), h=NULL",
    {
        n <- 10
        res <- c()
        for( i in 1:1000){
            pvals <- runif(n)
            statValue <- BJStat(pvals, l,h)
            res[i] <- BJExactPvalue(statValue, n, l, h)
        }
        expect_true(abs(mean(res)-0.5) < 0.01)
    }
)

test_that(
    "BJ statistics and critical, l=c(0,1), h=NULL",
    {
        n <- 10
        alpha <- 0.05
        critical <- BJExactCritical(n, alpha, l, h)
        res <- c()
        for( i in 1:10000){
            pvals <- runif(n)
            statValue <- BJStat(pvals, l,h)
            res[i] <- statValue<critical
        }
        expect_true(abs(mean(res)-0.05) < 0.01)
    }
)


test_that(
    "BJ test function, l=c(0,1), h=NULL",
    {
        ## BJ test
        BJplus <- BJTestGenerator(l=l,h=h)
        n <- 10
        res <- c()
        for( i in 1:10000){
            pvals <- runif(n)
            res[i] <- BJplus(pvals, alpha)
        }
        expect_true(abs(mean(res)-0.05) < 0.01)
    }
)

