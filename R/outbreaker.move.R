#' Movements of augmented data and parameters for outbreaker2
#'
#' This function moves all parameters for outbreaker2.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @param moves a set of movement functions stored in a named list as returned by \code{\link{outbreaker.create.moves}}
#' @param data a list of data items as returned by \code{outbreaker.data}
#' @param param a list of parameters as returned by \code{outbreaker.mcmc.init}
#' @param config a list of settings as returned by \code{outbreaker.config}
#' @param rand  a list of items as returned by \code{outbreaker.rand.vec}
#'
#' @return a potentially modified list of parameters as returned by \code{outbreaker.mcmc.init}
#' 
outbreaker.move <- function(moves, data, config, param, rand){
    ## get number of moves ##
    J <- length(moves)
    
    ## RUN MCMC ##
    for(i in seq.int(2, config$n.iter, 1)){
        for(j in seq.int(J)){
            param <- moves[[j]](data=data, param=param, config=config, rand=rand)
        }

        ## store outputs if needed
        if((i %% config$sample.every) == 0){
            param <- outbreaker.mcmc.store(param=param, data=data)
        }

    } # end of the chain

    ## RETURN ##
    return(param)
} # end outbreaker.move
