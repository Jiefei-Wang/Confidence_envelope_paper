test_that(
    "BJ statistics",
    {
        ## Test the power of the BJ statistics
        n <- 100
        alpha <- 0.5
        BJ <- BJGlobal()
        res <- c()
        for(i in 1:10000){
            pvalues <- runif(n)
            BJ$reset()
            BJ$addSamples(pvalues)
            res[i] <- BJ$hypothesisTest(alpha)
        }
        expect_true(abs(mean(res) - alpha) < 0.01)
    }
)


test_that(
    "BJ critical",
    {
        ## Test the critical function
        n <- 10
        alpha <- 0.055
        BJ <- BJGlobal()
        critical <- BJ$criticalFunc(alpha,n)
        res <- rep(0,100000)
        for(i in 1:100000){
            pvalues <- runif(n)
            BJ$reset()
            stat <- BJ$addSamples(pvalues)
            res[i] <- stat>=critical
        }
        expect_true(abs(mean(res) - alpha) < 0.005)
    }
)