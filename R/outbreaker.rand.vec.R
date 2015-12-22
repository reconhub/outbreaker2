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



#' Fast random number generation
#'
#' These functions use closure programming for fast generation of (logged) random numbers from various distributions.
#' \itemize{
#' \item \code{fast.log.runif} creates a function which generates 'n' values from Unif(0,1)
#' \item \code{fast.runif1} creates an optimized function generating a single value from Unif(0,1)
#' \item \code{fast.log.rnorm} creates a function which generates 'n' values from Norm(mu, sigma)
#' \item \code{fast.log.rnorm1} creates an optimized function generating a single value from Norm(mu, sigma)
#' }
#'
#' @param size.batch the size of the pre-generated vectors of values; larger batches lead to faster computations but require more RAM.
#'
#' @rdname fast.rand
#' @aliases fast.runif fast.log.runif fast.log.runif1 fast.log.rnorm fast.log.rnorm1
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @export
#'
#' @examples
#' x <- fast.log.runif(10)
#' environment(x)$values
#' x(6)
#' environment(x)$counter
#' x(4)
#' environment(x)$counter
#' x(3)
#' environment(x)$counter
#' environment(x)$values
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



#' @rdname fast.rand
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




## fast.log.rnorm <- function(mean=0, sd=1, size.batch=5e4){
##     ## initialize array
##     values <- log(rnorm(size.batch, mean=mean, sd=sd))

##     ## initialize counter
##     counter <- 0

##     ## create returned function
##     out <- function(n){
##         ## read from array if enough values
##         if((counter+n) <= size.batch){
##             counter <<- counter+n
##             return(values[seq.int(counter-n+1, counter)])
##         } else {
##             ## else, regenerate vector of values
##             if(n>size.batch) {
##                 size.batch <<- n
##             }
##             values <<- log(rnorm(size.batch, mean=mean, sd=sd))
##             counter <<- 0

##             return(out(n))
##         }
##     }

##     ## return output function
##     return(out)
## } # end outbreaker.log.rnorm



## #' @rdname fast.rand
## #' @export
## fast.log.rnorm1 <- function(mean=0, sd=1, size.batch=5e4){
##     ## initialize array
##     values <- log(rnorm(size.batch, mean=mean, sd=sd))

##     ## initialize counter
##     counter <- 0

##     ## create returned function
##     out <- function(mean=0, sd=1){
##         ## read from array if enough values
##         if(counter < size.batch){
##             counter <<- counter+1
##             return(values[counter])
##         } else {
##             ## else, regenerate vector of values
##             values <<- log(rnorm(size.batch, mean=mean, sd=sd))
##             counter <<- 0
##             return(out())
##         }
##     }

##     ## return output function
##     return(out)
## } # end fast.log.rnorm1
