## Fisher statistics
fisherStat <- function(x){
    x[x==0] <- 1e-10
    x[x==1] <- 1-1e-10
    -2*sum(log(x))
}
updataFisherStat <- function(statValue, samples, x){
    x[x==0] <- 1e-10
    x[x==1] <- 1-1e-10
    statValue - 2*sum(log(x))
}

fisherCritical <- function(alpha, n){
    qchisq(1-alpha, df = 2*n, lower.tail = TRUE)
}

FisherGlobal <- function(){
    .GlobalNullStat(statFunc = fisherStat, 
                        statUpdataFunc = updataFisherStat,
                        criticalFunc = fisherCritical)
}

## Stoffer statistics
stoufferStat <- function(x){
    x[x==0] <- 1e-10
    x[x==1] <- 1-1e-10
    -sum(qnorm(x))
}
updateStofferStat <- function(statValue, samples, x){
    x[x==0] <- 1e-10
    x[x==1] <- 1-1e-10
    statValue - qnorm(x)
}
stoufferCritical <- function(alpha,n){
    qnorm(1-alpha)*sqrt(n)
}

StofferGlobal <- function(){
    .GlobalNullStat(statFunc = stoufferStat, 
                        statUpdataFunc = updateStofferStat,
                        criticalFunc = stoufferCritical)
}


