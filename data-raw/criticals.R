
devtools::load_all()
library(foreach)
library(doParallel)
library(parallel)
cl <- makeCluster(10)
registerDoParallel(cl)
parallel::clusterEvalQ(cl, devtools::load_all())


# BJCriticalSpace <- new.env(parent = emptyenv())

for(n in 101:1000){
    alpha <- 0
    key <- getBJKey(n,alpha,TRUE)
    assign(key,0,envir = BJCriticalSpace)
    alpha <- 1
    key <- getBJKey(n,alpha,TRUE)
    assign(key,1,envir = BJCriticalSpace)

    alpha_list <- seq(0.01,0.99,0.01)
    final <- foreach(alpha=seq(0.01,0.99,0.01)) %dopar% {
        res <- genericCritical(
        pvalueFunc= BJExactPvalue,searchRange=c(0,1),
        n=n,alpha=alpha,
        oneSide = TRUE
        )
        res
    }
    for(i in seq_along(alpha_list)){
        alpha <- alpha_list[i]
        key <- getBJKey(n,alpha,TRUE)
        assign(key,final[[i]],envir = BJCriticalSpace)
    }
    message(n)
}

usethis::use_data(BJCriticalSpace, internal = TRUE, overwrite = TRUE)
