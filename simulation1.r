devtools::load_all()
source("setup.r",local =TRUE)
library(BiocParallel)
library(cp4p)
param <- SnowParam(workers = 20, type = "SOCK", progressbar =TRUE, RNGseed = 1, tasks=100)



n <- 1000
pi0 <- 0.5
n_null <- n*pi0
mu <- 2
trueNulls <- seq_len(n_null)
lambda_list <- seq(0,0.9,0.1)



nsim <- 1000

simulation1 <- function(x, lambda_list, pi_generator, generatePvals, n, n_null, mu, trueNulls){
    devtools::load_all(quiet = TRUE)
    pvals <- generatePvals(n, n_null, mu)
    pis <- lapply(lambda_list, pi_generator, pvals=pvals)
    story_pis <- lapply(lambda_list, story_generator)
    res_GCB <- c()
    res_true_null <- c()
  
    for(i in seq_along(pis)){
        res_true_null[i] <- pis[[i]](trueNulls)
        res_GCB[i] <- pi_BJ_fast(pis[[i]], pvals, 0.05, c(0,1), c(0,1))$upper
    }

    list(res_GCB, res_true_null)
}

res <- bplapply(1:nsim, simulation1, BPPARAM = param, 
lambda_list=lambda_list, 
pi_generator = pi_generator, 
generatePvals = generatePvals, n = n, n_null = n_null, mu = mu,
trueNulls = trueNulls)

## results are list of lists, 
## combine elelments of each list by row
all_res <- Reduce(bindList, res)

res_GCB <- all_res[[1]]
res_true_null <- all_res[[2]]

cover_prob <- colMeans(res_true_null<=res_GCB)
cover_prob

## ggplot to show the coverage probability as a function of lambda
library(ggplot2)
df <- data.frame(lambda=lambda_list, cover_prob=cover_prob)
ggplot(df, aes(x=lambda, y=cover_prob)) + geom_point() + geom_line() + theme_minimal() + labs(x="Lambda", y="Coverage Probability") + ylim(0.94,1) + geom_hline(yintercept=0.95, linetype="dashed", color = "black") + 
scale_y_continuous(breaks = seq(0.94, 1, by = 0.01))
## save
ggsave("coverage_prob.png")
