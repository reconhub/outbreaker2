#' Movements of augmented data and parameters for outbreaker2
#'
#' This function moves all parameters for outbreaker2.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @param param a list of parameters as returned by \code{outbreaker.create.mcmc}
#' @param rand  a list of items as returned by \code{outbreaker.rand.vec}
#' @param loglike  a list of log-likelihood functions with enclosed data as returned by \code{create.loglike}
#' @param priors  a list of prior functions with enclosed parameters as returned by \code{create.priors}
#' @param posteriors a list of posterior functions with enclosed likelihood and prior functions as
#' returned by \code{create.posteriors}
#'
#' @return a potentially modified list of parameters as returned by \code{outbreaker.create.mcmc}
#'
outbreaker.move <- function(moves, data, param, config, densities) {
    ## get number of moves ##
    J <- length(moves)

    ## RUN MCMC ##
    for (i in seq.int(2, config$n.iter, 1)) {
        ## move parameters / augmented data
        for (j in seq_len(J)) {
            ## safemode
            if (config$paranoid) {
                old.param <- param
            }

            ## move parameters
            param <- moves[[j]](param=param)

            ## safemode
            if (config$paranoid) {
                diagnostic <- look.for.trouble(param, data)
                if (!diagnostic$pass) {
                    stop(paste0(
                        "\n\n PARANOID MODE DETECTED AN ERROR WHILE FINDING IMPORTS:\n",
                        sprintf("at iteration %d, movement %d (%s) with the following diagnostics:\n%s\n",
                                i, j, names(moves)[j], diagnostic$msg)))
                }
            }
        }

        ## store outputs if needed
        if ((i %% config$sample.every) == 0) {
            param <- outbreaker.mcmc.store(param=param, densities=densities, step=i)
        }

    } # end of the chain

    ## RETURN ##
    return(param)
}
