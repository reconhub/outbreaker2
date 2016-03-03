#' Movements of augmented data and parameters for outbreaker2
#'
#' This function moves all parameters for outbreaker2.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @param data a list of data items as returned by \code{outbreaker.data}
#' @param param a list of parameters as returned by \code{outbreaker.mcmc.init}
#' @param config a list of settings as returned by \code{outbreaker.config}
#' @param rand  a list of items as returned by \code{outbreaker.rand.vec}
#'
#' @return a potentially modified list of parameters as returned by \code{outbreaker.mcmc.init}
#' 
outbreaker.move <- function(data, config, param, rand){
    ## RUN MCMC ##
    for(i in seq.int(2, config$n.iter, 1)){
        ## move infection dates ##
        if(config$move.t.inf){
            param <- move.t.inf(data=data, param=param, rand=rand)
        }

        ## move ancestries ##
        if(config$move.ances){
            param <- move.ances(data=data, param=param, config=config, rand=rand)
        }

        ## swap ancestries ##
        if(config$move.ances && config$move.t.inf){
            param <- move.swap.ances(data=data, param=param, config=config, rand=rand)
        }

        ## move mu ##
        if(config$move.mu){
            param <- move.mu(data=data, param=param, rand=rand)
        }

        ## store outputs if needed
        if((i %% config$sample.every) == 0){
            param <- outbreaker.mcmc.store(param=param, data=data)
        }

    } # end of the chain

    ## RETURN ##
    return(param)
} #end outbreaker.move
