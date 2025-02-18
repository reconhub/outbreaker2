## #' Movements of augmented data and parameters for outbreaker2
## #'
## #' This function moves all parameters for outbreaker2.
## #'
## #' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
## #'
## #'
## #' @param moves a list of movement functions as returned by \code{bind_moves}
## #'     (internal function)
## #'
## #' @param param a list of parameters as returned by
## #'     \code{create_param}
## #'
## #' @param densities a list containing lists of functions computing densities,
## #'     named: 'loglike' (log-likelihoods), 'priors' and 'posteriors'
## #'
## #' @return a potentially modified list of parameters as returned by
## #'     \code{create_param}
## #'
outbreaker_move <- function(moves, data, param_current,
                            param_store, config,
                            likelihoods, priors) {
  ## get number of moves ##
  J <- length(moves)

  ## Set up progress bar
  if(config$pb) {
    pb <- utils::txtProgressBar(min = 1, max = config$n_iter, style = 3)
  }

  ## RUN MCMC ##
  for (i in seq.int(2, config$n_iter, 1)) {
    ## move parameters / augmented data
    for (j in seq_len(J)) {
      ## move parameters
      param_current <- moves[[j]](param_current)

      ## safemode
      if (config$paranoid) {
        diagnostic <- look_for_trouble(param_current, param_store, data)
        if (!diagnostic$pass) {
          stop(paste0(
            "\n\n PARANOID MODE DETECTED AN ERROR WHILE FINDING IMPORTS:\n",
            sprintf("at iteration %d, movement %d (%s) with the following diagnostics:\n%s\n",
                    i, j, names(moves)[j], diagnostic$msg)))
        }
      }
    }

    ## store outputs if needed
    if ((i %% config$sample_every) == 0) {
      if(config$pb) {
        utils::setTxtProgressBar(pb, i)
      }
      param_store <- outbreaker_mcmc_store(param_current, param_store, data,
                                           config, likelihoods, priors, i)
    }

  } # end of the chain

  if(config$pb) {
    cat("\n")
  }

  ## output is a list of saved chain states
  return(param_store)
}
