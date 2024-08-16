
## fast algorithm
generatePvals <- function(n, n_null, a, b){
    pval_null <- runif(n_null)
    pval_alt <- rbeta(n - n_null, a,b)
    c(pval_null, pval_alt)
}
n <- 6
n_null <- 3
a <- 1
b <- 10
trueNulls <- seq_len(n_null)

alpha<-0.05
func <- pi_generator(0.5)

for(i in 1:1000){
    pvals <- generatePvals(n, n_null, a, b)
    res1 <- generalCE(pvals, func, alpha, BJTest)$upper
    res2 <- pi_BJ_fast(pvals, alpha, func)$upper
    stopifnot(all(res1 == res2))
}

func(1:5)
func(c(1,3,4,5,6))
