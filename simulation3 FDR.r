## This script is used to simulate the FDR and power of different methods under different pi0 values
## You MUST create your parallel backend beforehand, e.g.:
## param <- SnowParam(workers = 20)

# devtools::load_all()
library(generalExceedance)
library(cp4p)
library(BiocParallel)
source("setup.r",local =TRUE)
pi0_list <- seq(0,1,0.1)

FDP_table <- c()
conditional_FDP_table <- c()
power_table <- c()


for(pi0 in pi0_list){
    message("pi0: ", pi0)
    n <- 1000
    n_null <- n*pi0
    mu <- 2
    trueNulls <- seq_len(n_null)
    lambda <- 0.5
    alpha <- 0.05
    nsim <- 1000

    set.seed(1)
    pvalsList <- lapply(1:nsim, function(x) generatePvals(n, n_null, mu))


    simulation2 <- function(x, lambda, pi_generator, story_generator, generatePvals, n, n_null, mu, trueNulls, alpha, pvalsList){
        pvals <- pvalsList[[x]]
        pi_fun <- pi_generator(lambda, pvals)
        res_others <- c()
        res_GCB <- c()
    
        res_GCB <- generalExceedance::pi_BJ_fast(pi_fun, pvals, 0.05, c(0,1), c(0,1))$upper
        methods <- c("st.boot", "st.spline", "langaas", "jiang", "histo", "pounds", "abh","slim")
        # methods <- c("abh", "st.spline", "st.boot", "langaas")
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
            null_idx <- seq(n_null)
            alt_idx <- seq(n_null+1, n)
        }

        FDP <- c()
        power <- c()
        n_rejects <- c()
        for(i in seq_along(pi0_estimates)){
            # cur_alpha <- alpha_adj[i]
            cur_pi0 <- pi0_estimates[i]
            pvals_adj <- cp4p::adjust.p(pvals, pi0 = cur_pi0)$adjp$adjusted.p
            if(n_null == 0){
                FDP[i] <- 0
            }else{
                FDP[i] <- sum(pvals_adj[null_idx] <= alpha)/ max(1, sum(pvals_adj <= alpha))
            }
            if(n_null == n){
                power[i] <- 1
            }else{
                power[i] <- sum(pvals_adj[alt_idx] <= alpha)/ max(1, n-n_null)
            }
            n_rejects[i] <- sum(pvals_adj <= alpha)
        }

        names(FDP) <- names(power) <- names(pi0_estimates)

        list(FDP, power, pvals, n_rejects)
    }

    res <- bplapply(1:nsim, simulation2, BPPARAM = param, 
    lambda=lambda, 
    pi_generator = pi_generator, story_generator = story_generator, 
    generatePvals = generatePvals, n = n, n_null = n_null, mu = mu, trueNulls = trueNulls, alpha = alpha, pvalsList = pvalsList)

    ## results are list of lists, 
    ## combine elelments of each list by row
    all_res <- Reduce(bindList, res)

    res_FDP <- all_res[[1]]
    res_power <- all_res[[2]]
    res_pvals <- all_res[[3]]
    res_n_jects <- all_res[[4]]


    has_rejections <- res_n_jects > 0

    ave_FDP <- colMeans(res_FDP)
    ave_conditional_FDP <- colSums(res_FDP)/colSums(has_rejections)
    ave_power <- colMeans(res_power)

    # idx <- which(res_FDP[,1]>0)
    # pvals <- res_pvals[45,]

    stopifnot(!any(is.na(ave_FDP)))
    FDP_table <- cbind(FDP_table, ave_FDP)
    power_table <- cbind(power_table, ave_power)
}

colnames(FDP_table) <- paste0("pi0: ", pi0_list)
colnames(power_table) <- paste0("pi0: ", pi0_list)

save(FDP_table, power_table, file = "simulation3.RData")

## load data
load("simulation3.RData")

## save excel
library(openxlsx)
FDP_table2 <- as.data.frame(cbind(rownames(FDP_table), signif(FDP_table,2)))
power_table2 <- as.data.frame(cbind(rownames(power_table), signif(power_table,2)))
write.xlsx(FDP_table2, "FDP_table.xlsx")
write.xlsx(power_table2, "power_table.xlsx")

library(tidyverse)
FDP_plot <- t(FDP_table)|>
as.data.frame()|>
mutate(pi0 = pi0_list) |>
pivot_longer(-pi0, names_to = "method", values_to = "FDP") 

power_plot <- t(power_table)|>
as.data.frame()|>
mutate(pi0 = pi0_list) |>
pivot_longer(-pi0, names_to = "method", values_to = "power")

selections <- c("abh", "st.spline", "st.boot", "langaas", "jiang", "GCB")
## Make line plot to show FDP for each method as a function of pi0
library(ggplot2)
library(gridExtra)

FDP_plot|>
filter(method %in% selections) |>
mutate(method = factor(method, levels = selections))|>
filter(pi0>0)|>
ggplot(aes(x=pi0, y=FDP, group=method,shape = method)) + geom_point() + geom_line() + theme_minimal() + labs(x="pi0", y="FDR") + ylim(0,0.1) +
theme(axis.text.x = element_text(angle = 45, hjust = 1))+
geom_hline(yintercept=0.05, linetype="dashed", color = "black") +
scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5),
    labels = scales::percent_format(accuracy = 1))+ 
    theme(legend.position="none")

## save
ggsave("FDP_plot.png", width = 4, height = 4)


power_plot|>
filter(method %in% selections) |>
mutate(method = factor(method, levels = selections))|>
filter(pi0<1) |>
ggplot(aes(x=pi0, y=power, group=method,shape = method)) + geom_point() + geom_line() + theme_minimal() + labs(x="pi0", y="Power") + 
theme(axis.text.x = element_text(angle = 45, hjust = 1))+
scale_y_continuous(
    breaks = scales::pretty_breaks(n = 10),
    labels = scales::percent_format(accuracy = 1))

## save
ggsave("power_plot.png", width = 5, height = 4)






g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p1)

p3 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                         p2 + theme(legend.position="none"),
                         nrow=1),
             mylegend, nrow=2,heights=c(10, 1))