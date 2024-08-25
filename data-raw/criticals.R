library(BiocParallel)
# param <- SnowParam(workers = 20)
# BJCriticalSpace <- new.env(parent = emptyenv())
bptasks(param) <- 10000



n_list <- 1:5000
alpha <- 0.05
l <- c(0,1)
h <- c(0,1)
# final <- foreach(n=n_list) %dopar% {
#   BJExactCritical(n,alpha, l, h)
# }

res <- bplapply(n_list, function(n) {
  generalExceedance:::BJExactCritical(n,alpha, l, h)
}, BPPARAM = param)

devtools::load_all()
for(i in seq_along(n_list)){
    n <- n_list[i]
    key <- getBJKey(n,alpha, l, h)
    assign(key,res[[i]],envir = BJCriticalSpace)
}

usethis::use_data(BJCriticalSpace, internal = TRUE, overwrite = TRUE)
