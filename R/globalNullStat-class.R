
## GlobalNullStat class
## statFunc: function to calculate the statistic, large value means significance
## statUpdataFunc: function to update the statistic
## criticalFunc: function to calculate the critical value
.GlobalNullStat <- setRefClass(
    "GlobalNullStat",
    fields = list(
        statFunc = "ANY",
        statUpdataFunc = "ANY",
        criticalFunc = "ANY",
        samples = "numeric",
        statValue = "numeric"
    ),
    methods = list(
        initialize=function(statFunc=NULL, statUpdataFunc=NULL, criticalFunc){
            .self$statFunc <- statFunc
            .self$statUpdataFunc <- statUpdataFunc
            .self$criticalFunc <- criticalFunc
            .self$samples <- numeric(0)
            .self$statValue <- numeric(0)
        },
        addSamples = function(x){
            if (is.null(statUpdataFunc)) {
                ## add x to samples and update statValue
                samples <<- append(.self$samples, x)
                statValue <<- statFunc(samples)
            } else {
                if (is.null(statFunc)){
                    statValue <<- statUpdataFunc(statValue, samples, x)
                }else{
                    statValue <<- statFunc(x)
                }
                samples <<- append(.self$samples, x)
            }
            statValue
        },
        hypothesisTest = function(alpha){
            if(length(samples)==0)
                stop("No samples added!")
            statValue > criticalFunc(alpha, length(samples))
        },
        reset = function(){
            samples <<- numeric(0)
            statValue <<- numeric(0)
        }
    )
)
