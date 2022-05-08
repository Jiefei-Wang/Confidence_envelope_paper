BJStat <- function(x){
    exceedance::GKSStat(x, statName = "BJ+", pvalue = FALSE)$statValue
}
BJCritical <- function(alpha, n){
    exceedance::GKSCritical(n,alpha, statName = "BJ+")
}

BJGlobal <- function(){
    GlobalNullStat(statFunc = BJStat, criticalFunc = BJCritical)
}


StoufferStat <- function(pvalues){
    sum(qnorm(pvalues))
}

StoufferCritical <- function(n, alpha){
    qnorm(alpha)*sqrt(n)
}

StoufferGlobal <- function(){
    GlobalNullStat(statFunc = StoufferStat, criticalFunc = StoufferCritical)
}




