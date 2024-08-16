
devtools::load_all()
library(foreach)
library(doParallel)
library(parallel)
cl <- makeCluster(10)
registerDoParallel(cl)
parallel::clusterEvalQ(cl, devtools::load_all())


# BJCriticalSpace <- new.env(parent = emptyenv())

n_list <- 1:1000
alpha <- 0.05
l <- c(0,1)
h <- c(0,1)
final <- foreach(n=n_list) %dopar% {
  BJExactCritical(n,alpha, l, h)
}
for(i in seq_along(n_list)){
    n <- n_list[i]
    key <- getBJKey(n,alpha, l, h)
    assign(key,final[[i]],envir = BJCriticalSpace)
}

usethis::use_data(BJCriticalSpace, internal = TRUE, overwrite = TRUE)
