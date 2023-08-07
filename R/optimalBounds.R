findOptimalBounds <- function(globalTest, pvalues, groups, alpha, FC = 2){
    nGroups <- length(groups)
    groupSize <- length(groups[[1]])
    potential <- rep(0, nGroups)
    ## Find the number of hypotheses that are significant in each group
    for(i in seq_len(nGroups)){
        cutoff <- rep(0, nGroups)
        for(j in seq_len(groupSize)){
            cutoff[i] <- j
            if(FalseCoverage(globalTest, pvalues, groups, cutoff, alpha) == 0){
                potential[i] <- j
                break
            }
        }
    }
    cutoff <- rep(0, nGroups)
    lastFC <- 0
    tempFC <- rep(0, nGroups)
    blockList <- c()
    while(TRUE){
        candidateIdx <- c()
        for(i in seq_len(nGroups)){
            curCutoff <- cutoff
            if(curCutoff[i]==groupSize)
                next
            curCutoff[i] <- curCutoff[i] + 1
            curFC <- FalseCoverage(globalTest, pvalues, groups, curCutoff, alpha)
            tempFC[i] <- curFC
            if(curFC <= FC){
                candidateIdx <- c(candidateIdx, i)
            }
        }
        ## Exclude the indices that are already in the block list
        candidateIdx <- setdiff(candidateIdx, blockList)
        if(length(candidateIdx) == 0){
            break
        }else{
            indices <- candidateIdx
            idx <- which.max(potential[indices] - cutoff[indices])
            idx <- indices[idx]
            cutoff[idx] <- cutoff[idx] + 1
            lastFC <- tempFC[idx]

            ## If incrementing a cutoff does not change the false coverage, add it to the block list
            curCutoff <- cutoff
            curCutoff[idx] <- groupSize
            curFC <- FalseCoverage(globalTest, pvalues, groups, curCutoff, alpha)
            if(curFC==lastFC){
                blockList <- c(blockList, idx)
            }
        }
    }
    list(cutoff=cutoff, FC=lastFC)
}
