################
## selectChains
################
selectChains <- function(x, select="visual", alpha=0.001, ...){
    ## CHECKS ##
    if(!is.list(x)) stop("x should be a list as output by outbreaker / outbreaker.parallel")
    if(x$n.runs==1) { # return x unchanged if only one run
        return(x)
    }
    if(is.character(select)) select <- match.arg(select, c("auto","visual"))


    ## SELECTION BASED ON VISUAL INSPECTION OF THE CHAINS ##
    if(is.character(select) && select=="visual"){
        repeat{
            plotChains(x, ...)
            if(x$n.runs==1) {
                cat("Only one run left - exiting.\n")
                break
            }
            cat("Indicate the runs to remove from the results (0 for exit): ")
            answer <- suppressWarnings(as.integer(readLines(n = 1)))
            if(is.na(answer) || answer==0 || answer>x$n.runs) break
            x <- selectChains(x, select=setdiff(1:x$n.runs, answer))
            cat(paste("(removed run ", answer,")\n",sep=""))
        }
    }


    ## SELECTION BASED ON AUTOMATIC PROCEDURE ##
    ## idea:
    ## 1) make ANOVA
    ## if significant:
    ##   2a) remove the run with smallest log-post
    ##   2b) go to 1)
    ## if not: exit
    if(is.character(select) && select=="auto"){
        repeat{
            if(x$n.runs==1) {
                cat("Only one run left - exiting.\n")
                break
            }
            pval <- anova(lm(x$chains$post ~ factor(x$chains$run)))$"Pr(>F)"[1]
            if(pval >= alpha) break
            toRemove <- which.min(tapply(x$chains$post, factor(x$chains$run), mean))
            x <- selectChains(x, select=setdiff(1:x$n.runs, toRemove))
            cat(paste("(removed run ", toRemove,")\n",sep=""))
        }
    }


    ## SELECTION PROVIDED AS NUMBERS ##
    if(is.numeric(select) || is.integer(select)){
        x$chains <- x$chains[x$chains$run %in% select, , drop=FALSE]
        x$chains$run <- as.integer(factor(x$chains$run))
        x$n.runs <- length(select)
    }

    return(x)

} # end selectChains
