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
outbreaker.move <- function(moves, data, param, config, densities, rand){
    ## get number of moves ##
    J <- length(moves)

    ## RUN MCMC ##
    for(i in seq.int(2, config$n.iter, 1)){
        ## move parameters / augmented data
        for(j in seq.int(J)){
            ## safemode
            if(config$paranoid){
                old.param <- param
            }

            ## move parameters
            param <- moves[[j]](param=param, config=config, densities=densities, rand=rand)

            ## safemode
            if(config$paranoid){
                diagnostic <- look.for.trouble(param, data)
                if(!diagnostic$pass){
                    stop(paste("\n\n PARANOID MODE DETECTED AN ERROR:\n",
                               "at iteration ", i, ", ",
                               "movement ", j, " (", names(moves)[j], ")",
                               " with the following diagnostics:\n", diagnostic$msg, "\n\n",
                               sep=""))
                }
            }
        }

        ## store outputs if needed
        if((i %% config$sample.every) == 0){
            param <- outbreaker.mcmc.store(param=param, data=data, densities=densities, step=i)
        }

    } # end of the chain

    ## RETURN ##
    return(param)
} # end outbreaker.move
