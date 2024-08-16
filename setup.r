

generatePvals <- function(n, n_null, mu){
  pval_null <- runif(n_null)
  pval_alt <- 1-pnorm(rnorm(n - n_null, mu))
  c(pval_null, pval_alt)
}


pi_generator <- function(lambda, pvals){
  force(pvals)
  function(idx){
    sum(pvals[idx] >lambda)/(length(pvals)*(1-lambda))
  }
}


story_generator <- function(lambda){
  function(pvals){
    sum(pvals >lambda)/(length(pvals)*(1-lambda))
  }
}



bindList <- function(x,y){
  ## rbind list elements
  for(i in seq_along(x)){
    x[[i]] <- rbind(x[[i]], y[[i]])
  }
  x
}

MSE <- function(res, pi0){
  colMeans((res-pi0)^2)
}


my_round = function(x, n=2) {
  max(abs(round(x, n)), abs(signif(x, 1))) * sign(x)
  }
