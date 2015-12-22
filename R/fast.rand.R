#' Fast random number generation
#'
#' These functions use closure programming for fast generation of random numbers from various distributions.
#' \itemize{
#' \item \code{fast.log.runif} creates a function which generates 'n' (logged) values from Unif(0,1)
#' \item \code{fast.runif1} creates an optimized function generating a single (logged) value from Unif(0,1)
#' \item \code{fast.rnorm} creates a function which generates 'n' values from Norm(mean, sd)
#' \item \code{fast.rnorm1} creates an optimized function generating a single value from Norm(mean, sd)
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
} # end fast.log.runif





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





#' @rdname fast.rand
#' @export
fast.rnorm <- function(mean=0, sd=1, size.batch=5e4){
    ## initialize array
    values <- rnorm(n=size.batch, mean=mean, sd=sd)

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
            values <<- rnorm(size.batch, mean=mean, sd=sd)
            counter <<- 0

            return(out(n))
        }
    }

    ## return output function
    return(out)
} # end fast.rnorm





#' @rdname fast.rand
#' @export
fast.rnorm1 <- function(mean=0, sd=1, size.batch=5e4){
    ## initialize array
    values <- rnorm(n=size.batch, mean=mean, sd=sd)

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
            values <<- rnorm(n=size.batch, mean=mean, sd=sd)
            counter <<- 0
            return(out())
        }
    }

    ## return output function
    return(out)
} # end fast.rnorm1
