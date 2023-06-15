FalseCoverage <- function(globalTest,pvalues,groups,cutoff, alpha){
    ## Remove the groups with cutoff == 0
    groups_filtered <- groups[cutoff != 0]
    cutoff_filtered <- cutoff[cutoff != 0]

    ## change the group indices to the indices in sorted p-value
    groups_filtered <- lapply(groups_filtered, function(x)
        rank(-pvalues)[x])


    groupsTopPIdx <- list()
    for(i in seq_along(groups_filtered)){
        l <- length(groups_filtered[[i]]) - cutoff_filtered[i] + 1
        group_idx <- sort(groups_filtered[[i]])
        ## Get the largest l p-value in each group
        groupsTopPIdx[[i]] <- group_idx[seq_len(l)]
    }

    n<- length(pvalues)
    k <- length(groups_filtered)
    p_sort <- sort(pvalues, decreasing = TRUE)
    ## Number of false coverage in groups
    for(FC in rev(seq_len(k))){
        for (j in seq_len(choose(k, FC))) {
            ## Generate the next combination
            if (j == 1) {
                cbn <- seq_len(FC)
            } else {
                cbn <- gen.next.cbn(cbn, k)
            }
            ## Get the p-values in the combination
            p1_idx <- unlist(groupsTopPIdx[cbn])
            p1 <- p_sort[p1_idx]

            ## The pvlaues not in the combination
            p_rest <- p_sort[-p1_idx]

            ## Initialize samples
            globalTest$reset()
            globalTest$addSamples(p1)
            ## Append the p-values not in the combination
            for (q in 0:(n-length(p1))){
                if(q!=0){
                    globalTest$addSamples(p_rest[q])
                }
                ## If test is not significant, stop and return the combination
                if(!globalTest$hypothesisTest(alpha)){
                    #c(p1_idx, seq_len(n)[-p1_idx][1:q]))
                    return(FC)
                }
            }
        }
    }
    0
}
