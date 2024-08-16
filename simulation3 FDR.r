devtools::load_all()
source("setup.r",local =TRUE)
library(BiocParallel)
library(cp4p)
param <- SnowParam(workers = 20, type = "SOCK", progressbar =TRUE, RNGseed = 1, tasks=100)

pi0_list <- seq(0,1,0.2)

FDP_table <- c()
power_table <- c()


for(pi0 in pi0_list){
    message("pi0: ", pi0)
    n <- 1000
    n_null <- n*pi0
    mu <- 2
    trueNulls <- seq_len(n_null)
    lambda <- 0.5


    nsim <- 1000
    simulation2 <- function(x, lambda, pi_generator, story_generator, generatePvals, n, n_null, mu, trueNulls){
        devtools::load_all(quiet = TRUE)
        pvals <- generatePvals(n, n_null, mu)
        pi_fun <- pi_generator(lambda, pvals)
        res_others <- c()
        res_GCB <- c()
    
        res_GCB <- pi_BJ_fast(pi_fun, pvals, 0.05, c(0,1), c(0,1))$upper
        methods <- c("st.boot", "st.spline", "langaas", "jiang", "histo", "pounds", "abh","slim")
        res_others<- sapply(methods, function(method) cp4p::estim.pi0(pvals, pi0.method = method)[[1]])

        ## all estimates, clamped to [alpha,1], 
        ## anything less than alpha means rejecting all
        pi0_estimates <- c(GCB = res_GCB, res_others)
        pi0_estimates[pi0_estimates<alpha] <- alpha
        pi0_estimates[pi0_estimates>1] <- 1

        # pvals_adj <- p.adjust(pvals, method = "BH")
        # alpha_adj <- alpha/pi0_estimates


        
        ## index of null and alternative
        if(n_null == 0){
            null_idx <- integer(0)
            alt_idx <- seq(n)
        }else if(n_null == n){
            null_idx <- seq(n)
            alt_idx <- integer(0)
        }else{
            null_idx <- sample(seq(n), n_null)
            alt_idx <- seq(n_null+1, n)
        }

        FDP <- c()
        power <- c()
        for(i in seq_along(pi0_estimates)){
            # cur_alpha <- alpha_adj[i]
            cur_pi0 <- pi0_estimates[i]
            pvals_adj <- adjust.p(pvals, pi0 = cur_pi0)$adjp$adjusted.p
            if(n_null == 0){
                FDP[i] <- 0
            }else{
                FDP[i] <- mean(pvals_adj[null_idx] <= alpha)/ max(1, sum(pvals <= alpha))
            }
            if(n_null == n){
                power[i] <- 1
            }else{
                power[i] <- mean(pvals_adj[alt_idx] <= alpha)
            }
        }

        names(FDP) <- names(power) <- names(pi0_estimates)

        list(FDP, power)
    }

    opt <- bpoptions(RNGseed = 1)
    res <- bplapply(1:nsim, simulation2, BPPARAM = param, 
    lambda=lambda, 
    pi_generator = pi_generator, story_generator = story_generator, 
    generatePvals = generatePvals, n = n, n_null = n_null, mu = mu, trueNulls = trueNulls, BPOPTIONS = opt)

    ## results are list of lists, 
    ## combine elelments of each list by row
    all_res <- Reduce(bindList, res)

    res_FDP <- all_res[[1]]
    res_power <- all_res[[2]]

    ave_FDP <- colMeans(res_FDP)
    ave_power <- colMeans(res_power)

    FDP_table <- cbind(FDP_table, ave_FDP)
    power_table <- cbind(power_table, ave_power)
}

colnames(FDP_table) <- paste0("pi0: ", pi0_list)
colnames(power_table) <- paste0("pi0: ", pi0_list)

## save excel
library(openxlsx)
write.xlsx(FDP_table, "FDP_table.xlsx")
write.xlsx(power_table, "power_table.xlsx")


