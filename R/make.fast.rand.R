#' Fast random number generation
#'
#' These functions use pre-generation of large random vectors to speed up random number generation.
#' \itemize{
#' \item \code{make.fast.rand} creates a function which generates 'n' values from a given distribution \code{f}.
#' \item \code{make.fast.runif1} creates an optimized function generating a single value from a given distribution \code{f}.
#' }
#'
#' @param ... arguments to be passed to a function implementing random number generation, e.g. \code{runif} or \code{rnorm}.
#' @param f a function implementing random number generation, e.g. \code{runif} or \code{rnorm}.
#' @param batch.size the size of the pre-generated vectors of values; larger batches lead to faster computations but require more RAM.
#' @param log a logical indicating whether generated values should be logged; defaults to FALSE.
#'
#' @rdname make.fast.rand
#' @aliases make.fast.rand make.fast.rand1
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}, Rich Fitzjohn
#'
#' @export
#'
#' @importFrom compiler cmpfun
#'
#' @examples
#' ## from Unif(0,1)
#' x <- make.fast.rand()
#' environment(x)$values
#' x(6)
#' environment(x)$counter
#' x(4)
#' environment(x)$counter
#' x(3)
#' environment(x)$counter
#' environment(x)$values
#'
#' ## check with microbenchmark
#' if(require(microbenchmark)){
#' x1 <- make.fast.rand1()
#' plot(microbenchmark(x1(), runif(1), times=1e3L), ylim=c(1000,3000))
#' }
#'
make.fast.rand <- function(..., f=runif, batch.size=5e4, log=FALSE){
    ## initialize array
    values <- f(batch.size, ...)
    if(log) values <- log(values)

    ## initialize counter
    counter <- 0

    ## create returned function
    out <- function(n){
        ## read from array if enough values
        if((counter+n) <= batch.size){
            counter <<- counter+n
            return(values[seq.int(counter-n+1, counter)])
        } else {
            ## else, regenerate vector of values
            if(n>batch.size) {
                batch.size <<- n
            }
            values <<- f(batch.size, ...)
            if(log) values <<- log(values)
            counter <<- 0

            return(out(n))
        }
    }

    ## return output function
    return(out)
} # end make.fast.rand

make.fast.rand <- compiler::cmpfun(make.fast.rand)



#' @rdname make.fast.rand1
#' @export
make.fast.rand1 <- function(..., f=runif, batch.size=5e4, log=FALSE){
    ## initialize array
    values <- f(batch.size, ...)
    if(log) values <- log(values)

    ## initialize counter
    counter <- 0

    ## create returned function
    out <- function(){
        ## read from array if enough values
        if(counter < batch.size){
            counter <<- counter+1
            return(values[counter])
        } else {
            ## else, regenerate vector of values
            values <<- f(batch.size, ...)
            if(log) values <<- log(values)
            counter <<- 0
            return(out())
        }
    }

    ## return output function
    return(out)
} # end make.fast.rand1

make.fast.rand1 <- compiler::cmpfun(make.fast.rand1)


