## This script is used to estimate the coverage probability
## You MUST create your parallel backend beforehand, e.g.:
## param <- SnowParam(workers = 20)
source("setup.r",local =TRUE)
library(BiocParallel)
library(cp4p)



n <- 1000
pi0 <- 0.5
n_null <- n*pi0
mu <- 2
trueNulls <- seq_len(n_null)
lambda_list <- seq(0,0.9,0.1)



nsim <- 1000
set.seed(1)
pvalsList <- lapply(1:nsim, function(x) generatePvals(n, n_null, mu))

simulation1 <- function(x, lambda_list, pi_generator, generatePvals, n, n_null, mu, trueNulls, pvalsList){
    pvals <- pvalsList[[x]]
    pis <- lapply(lambda_list, pi_generator, pvals=pvals)
    story_pis <- lapply(lambda_list, story_generator)
    res_GCB <- c()
    res_true_null <- c()
  
    for(i in seq_along(pis)){
        res_true_null[i] <- pis[[i]](trueNulls)
        res_GCB[i] <- generalExceedance::pi_BJ_fast(pis[[i]], pvals, 0.05, c(0,1), c(0,1))$upper
    }

    list(res_GCB, res_true_null)
}

res <- bplapply(1:nsim, simulation1, BPPARAM = param, 
lambda_list=lambda_list, 
pi_generator = pi_generator, 
generatePvals = generatePvals, n = n, n_null = n_null, mu = mu,
trueNulls = trueNulls, pvalsList=pvalsList)

## results are list of lists, 
## combine elelments of each list by row
all_res <- Reduce(bindList, res)

res_GCB <- all_res[[1]]
res_true_null <- all_res[[2]]

cover_prob <- colMeans(res_true_null<=res_GCB)


cover_prob <- as.data.frame(t(cover_prob))
colnames(cover_prob) <- paste0("lambda=",lambda_list)
cover_prob <- cbind("Cover Prob", cover_prob)

library(openxlsx)
write.xlsx(cover_prob, "cover_prob.xlsx")
