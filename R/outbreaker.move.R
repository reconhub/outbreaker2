#' Movements of augmented data and parameters for outbreaker2
#'
#' This function moves all parameters for outbreaker2.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @inheritParams outbreaker.create.mcmc
#'
#' @param moves a list of movement functions as returned by \code{create.moves} (internal function)
#'
#' @param param a list of parameters as returned by \code{outbreaker.create.mcmc}
#'
#' @param densities a list containing lists of functions computing densities, named: 'loglike' (log-likelihoods), 'priors' and 'posteriors'
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

    ## output is a list of saved chain states
    return(param)
}
