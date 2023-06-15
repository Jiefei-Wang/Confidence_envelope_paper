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



findNullSets <- function(pvalues, globalTest, alpha){
    hypos <- names(pvalues)
    if (is.null(hypos)) hypos <- seq_along(pvalues)
    results <- list()
    n <- length(hypos)
    for(k in rev(seq_along(hypos))){
        for (i in seq_len(choose(n, k))){
            if (i == 1)
                cbn <- 1:k
            else
                cbn <- gen.next.cbn(cbn, n)
            psubset <- pvalues[hypos[cbn]]
            globalTest$reset()
            globalTest$addSamples(psubset)
            if (!globalTest$hypothesisTest(alpha)) {
                results <- c(results, list(hypos[cbn]))
            }
        }
    }
    results
}

maxLoss <- function(lossFunc, nullSets){
    loss <- lossFunc(integer(0))
    for(i in nullSets){
        loss <- max(loss, lossFunc(i))
    }
    loss
}

minLoss <- function(lossFunc, nullSets){
    loss <- lossFunc(integer(0))
    ##loss <- .Machine$integer.max
    for(i in nullSets){
        loss <- min(loss, lossFunc(i))
    }
    loss
}


# GE <- function(pvalues, alpha, lossFunc, testFunc, nullSets = NULL){
#     if (is.null(nullSets))
#         nullSets <- findNullSets(pvalues, testFunc, alpha)
#     maxLoss(lossFunc, nullSets)
# }

GE <- function(pvalues, alpha, globalTest, lossFunc, nullSets = NULL){
    if (is.null(nullSets))
        nullSets <- findNullSets(pvalues, globalTest, alpha)
    list(
        min = minLoss(lossFunc, nullSets),
        max = maxLoss(lossFunc, nullSets)
    )
}
