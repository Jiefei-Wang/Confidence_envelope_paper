devtools::load_all()
source("setup.r",local =TRUE)
library(BiocParallel)
library(cp4p)
param <- SnowParam(workers = 20, type = "SOCK", progressbar =TRUE, RNGseed = 1, tasks=100)

pi0_list <- seq(0,1,0.2)

ave_table <- c()
SSE_table <- c()


for(pi0 in pi0_list){
    message("pi0: ", pi0)
    n <- 1000
    n_null <- n*pi0
    mu <- 2
    trueNulls <- seq_len(n_null)
    lambda_list <- c(0, 0.3, 0.5, 0.7, 0.9)


    nsim <- 1000
    simulation2 <- function(x, lambda_list, pi_generator, story_generator, generatePvals, n, n_null, mu, trueNulls){
        devtools::load_all(quiet = TRUE)
        pvals <- generatePvals(n, n_null, mu)
        pis <- lapply(lambda_list, pi_generator, pvals=pvals)
        story_pis <- lapply(lambda_list, story_generator)
        res_others <- c()
        res_true_null <- c()
        res_GCB <- c()
        res_story <- c()
    
        for(i in seq_along(pis)){
            res_story[i] <- story_pis[[i]](pvals)
            res_true_null[i] <- pis[[i]](trueNulls)
            res_GCB[i] <- pi_BJ_fast(pis[[i]], pvals, 0.05, c(0,1), c(0,1))$upper
        }
        methods <- c("st.boot", "st.spline", "langaas", "jiang", "histo", "pounds", "abh","slim")
        res_others<- sapply(methods, function(method) cp4p::estim.pi0(pvals, pi0.method = method)[[1]])

        list(res_true_null, res_GCB, res_others, res_story)
    }

    res <- bplapply(1:nsim, simulation2, BPPARAM = param, 
    lambda_list=lambda_list, 
    pi_generator = pi_generator, story_generator = story_generator, 
    generatePvals = generatePvals, n = n, n_null = n_null, mu = mu, trueNulls = trueNulls)

    ## results are list of lists, 
    ## combine elelments of each list by row
    all_res <- Reduce(bindList, res)

    res_true_null <- all_res[[1]]
    res_GCB <- all_res[[2]]
    res_others <- all_res[[3]]
    res_story <- all_res[[4]]

    ave_true_null <- colMeans(res_true_null)
    ave_GCB <- colMeans(res_GCB)
    ave_others <- colMeans(res_others)
    ave_story <- colMeans(res_story)

    names(ave_true_null) <- paste0("true null estimator ", lambda_list)
    names(ave_GCB) <- paste0("GCB estimator ", lambda_list)
    names(ave_story) <- paste0("story estimator ", lambda_list)


    ## SSE
    SSE_true_null <- SSE(res_true_null, pi0)
    SSE_GCB <- SSE(res_GCB, pi0)
    SSE_others <- SSE(res_others, pi0)
    SSE_story <- SSE(res_story, pi0)

    names(SSE_true_null) <- paste0("true null estimator ", lambda_list)
    names(SSE_GCB) <- paste0("GCB estimator ", lambda_list)
    names(SSE_story) <- paste0("story estimator ", lambda_list)

    ave_table <- cbind(ave_table, c(ave_true_null, ave_GCB, ave_others, ave_story))

    SSE_table <- cbind(SSE_table, c(SSE_true_null, SSE_GCB, SSE_others, SSE_story))
}

colnames(ave_table) <- paste0("pi0: ", pi0_list)
colnames(SSE_table) <- paste0("pi0: ", pi0_list)


## save excel
library(openxlsx)
write.xlsx(round(ave_table,3), "ave_table.xlsx")
write.xlsx(signif(SSE_table,2), "SSE_table.xlsx")

select <- c("st.boot", "st.spline", "langaas", "jiang",
"histo", "pounds", "abh", "slim", "GCB estimator 0",
"GCB estimator 0.3", "GCB estimator 0.5", "GCB estimator 0.7",
"GCB estimator 0.9"
)

SSE_table$row_ave <- apply(SSE_table, 1, mean)

write.xlsx(round(ave_table[select,],3), "ave_table_paper.xlsx")
write.xlsx(signif(SSE_table[select,],2), "SSE_table_paper.xlsx")

