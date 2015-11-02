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
outbreaker.rand.vec <- function(config){
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
