## Function definition
gen.next.cbn <- function(cbn, n){
    ## Generates the combination that follows the one provided as input
    cbn.bin      <- rep(0, n)
    cbn.bin[cbn] <- 1
    if (tail(cbn.bin, 1) == 0){
        ind <- tail(which(cbn.bin == 1), 1)
        cbn.bin[c(ind, ind+1)] <- c(0, 1)
    }else{
        ind <- 1 + tail(which(diff(cbn.bin) == -1), 1)
        nb  <- sum(cbn.bin[-c(1:ind)] == 1)
        cbn.bin[c(ind-1, (n-nb+1):n)] <- 0
        cbn.bin[ind:(ind+nb)]         <- 1
    }
    cbn <- which(cbn.bin == 1)
    cbn
}


findNullSets <- function(pvalues, testFunc, alpha){
    hypos <- names(pvalues)
    if (is.null(hypos)) hypos <- seq_along(pvalues)
    results <- list()
    for(k in length(hypos):1){
        for (i in seq_len(choose(n, k))){
            if (i == 1)
                cbn <- 1:k
            else
                cbn <- gen.next.cbn(cbn, n)
            if (testFunc(pvalues[hypos[cbn]]) > alpha) {
                results <- c(results, list(hypos[cbn]))
            }
        }
    }
    results
}

maxLoss <- function(lossFunc, nullSets){
    loss <- lossFunc(NULL)
    for(i in nullSets){
        loss <- max(loss, lossFunc(i))
    }
    loss
}

minLoss <- function(lossFunc, nullSets){
    if(length(nullSets)==0)
        return(0)
    loss <- .Machine$integer.max
    for(i in nullSets){
        loss <- min(loss, lossFunc(i))
    }
    loss
}


GE <- function(pvalues, alpha, lossFunc, testFunc, nullSets = NULL){
    if (is.null(nullSets))
        nullSets <- findNullSets(pvalues, testFunc, alpha)
    maxLoss(lossFunc, nullSets)
}


# I <- list(1:5,6:10,11:12)
# K <- c(2,2,1)
GCS_FP_param <- function(pvalues, I, globalStat) {
    p_sort <- sort(pvalues)
    p_rank <- rank(pvalues)
    I_rank <- rep(list(NULL), length(I))
    for(i in seq_along(I)){
        I_rank[[i]] <- sort(p_rank[I[[i]]])
    }
    criticals <- new.env(parent = emptyenv())

    list(
        p_sort = p_sort,
        p_rank = p_rank,
        I_rank = I_rank,
        criticals = criticals,
        globalStat = globalStat
    )
}

GCS_FP <- function(param, k, alpha){
    p_sort = param$p_sort
    p_rank = param$p_rank
    I_rank = param$I_rank
    criticals <- param$criticals
    globalStat <- param$globalStat
    m <- length(p_sort)

    k_zeros <- sum(k <= 0)
    exclude_k <- which(k <= 0 | k > lengths(I_rank))
    if(length(exclude_k)!=0){
        I_rank <- I_rank[-exclude_k]
        k <- k[-exclude_k]
    }
    if(length(k) == 0)
        return(k_zeros)

    I_len <- length(I_rank)
    I_in_sets <- rep(list(NULL), I_len)
    for(i in seq_along(I_rank)){
        I_in_sets[[i]] <- tail(I_rank[[i]], k[i])
    }

    for(n in rev(seq_len(I_len))){
        for (i in seq_len(choose(I_len, n))){
            if (i == 1)
                cbn <- 1:n
            else
                cbn <- gen.next.cbn(cbn, I_len)
            I_in <- unlist(I_in_sets[cbn])
            I_ex <- setdiff(seq_len(m), I_in)
            curStat <- stat(globalStat, p_sort[I_in])
            n_in <- length(I_in)
            for(j in rev(0:length(I_ex))){
                if(n_in != 0) {
                    critical_name <- as.character(n_in)
                    if(is.null(criticals[[critical_name]])){
                        criticals[[critical_name]] <-
                            critical(globalStat, alpha, n_in)
                    }
                    critical <- criticals[[critical_name]]

                    if(curStat > critical){
                        # browser()
                        return(n + k_zeros)
                    }
                }
                if(j!=0){
                    curStat <- updateStat(globalStat, curStat, p_sort[I_ex[j]])
                    n_in <- n_in + 1L
                }
            }
        }
    }
    return(k_zeros)
}


