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

    ## CREATE FAST RANDOM NUMBER GENERATORS ##
    out$log.runif1 <- make.fast.rand1(f=runif, config$batch.size, log=TRUE)
    out$mu.rnorm1 <- make.fast.rand1(f=rnorm, config$batch.size, mean=0, sd=config$sd.mu, log=FALSE)

    ## RETURN ##
    return(out)
} # end outbreaker.rand.vec
