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



#' Fast random uniform number generation
#'
#' These functions use closure programming for fast generation of (logged) random numbers from a uniform distribution on [0;1].
#' \itemize{
#' \item \code{fast.log.runif} creates a function which generates 'n' values
#' \item \code{fast.runif1} creates an optimized function generating a single value
#' }
#'
#' @param size.batch the size of the pre-generated vectors of values; larger batches lead to faster computations but require more RAM.
#'
#' @rdname fast.runif
#' @aliases fast.log.runif fast.log.runif1
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @export
#'
fast.log.runif <- function(size.batch=5e4){
    ## initialize array
    values <- log(runif(size.batch))

    ## initialize counter
    counter <- 0

    ## create returned function
    out <- function(n){
        ## read from array if enough values
        if((counter+n) <= size.batch){
            counter <<- counter+n
            return(values[seq.int(counter-n+1, counter)])
        } else {
            ## else, regenerate vector of values
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
} # end outbreaker.log.runif



#' @rdname fast.runif
#' @export
fast.log.runif1 <- function(size.batch=5e4){
    ## initialize array
    values <- log(runif(size.batch))

    ## initialize counter
    counter <- 0

    ## create returned function
    out <- function(){
        ## read from array if enough values
        if(counter < size.batch){
            counter <<- counter+1
            return(values[counter])
        } else {
            ## else, regenerate vector of values
            values <<- log(runif(size.batch))
            counter <<- 0
            return(out())
        }
    }

    ## return output function
    return(out)
} # end fast.log.runif1
