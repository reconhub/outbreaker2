#' Random number pre-generation for outbreaker2
#'
#' This function generates pre-drawn vectors of random numbers to be used in outbreaker2.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @rdname outbreaker.rand.vec
#'
#' @param config a list of settings as returned by \code{outbreaker.config}
#'
#' @export
#'
outbreaker.rand.vec <- function(config)
{
    ## CREATE OUTPUT ##
    out <- list()

    ## FILL IN ARRAYS OF RANDOM NUMBERS ##
    out$acc.t.inf <- log(runif(config$n.iter))
    out$acc.ances <- log(runif(config$n.iter))
    out$mu <- rnorm(config$n.iter, mean=0, sd=config$sd.mu)
    out$acc.mu <- log(runif(config$n.iter))

    ## RETURN ##
    return(out)
} # end outbreaker.rand.vec



counter <- function(batch.size=10){
    i <- 0
    out <- function(){
        if(i<batch.size){
            i <<- i+1
        } else {
            i <<- 0
        }
        return(i)
    }
    return(out)
}


outbreaker.log.runif <- function(size.batch=10){
    ## initialize array
    values <- log(runif(size.batch))

    ## initialize counter
    counter <- 0

    ## create returned function
    out <- function(n){
        ## read from array if enough values
        if((counter+n) <= size.batch){
            cat("\n reading values from vector\n")
            on.exit(return(values[seq.int(counter+1, counter+n)]))
            counter <<- counter+n
        } else {
            ## regenerate vector of values
            cat("\n generating vector\n")
            if(n>size.batch) {
                size.batch <<- n
            }
            values <<- log(runif(size.batch))
            counter <<- 0

            return(out(n))
        }
    }

    ## return output function
    return(out)
}
