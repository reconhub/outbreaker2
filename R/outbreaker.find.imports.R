#' Detection of imported cases for outbreaker2
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
outbreaker.find.imports <- function(moves, data, config, param, rand){
    ## return if not import ##
    if(!config$detect.import) return(list(config=config, param=param))


    ## store initial param values ##
    ini.param <- param

    ## get number of moves ##
    J <- length(moves)

    ## create matrix of individual influences ##
    n.measures <- floor(config$n.iter-1000/config$sample.every)
    influences <- matrix(0, ncol=data$N, nrow=n.measures)
    counter <- 1L

    ## RUN MCMC ##
    for(i in seq.int(2, config$n.iter, 1)){
        ## move parameters / augmented data
        for(j in seq.int(J)){
            ## safemode
            if(config$paranoid){
                old.param <- param
            }

            ## move parameters
            param <- moves[[j]](data=data, param=param, config=config, rand=rand)

            ## safemode
            if(config$paranoid){
                diagnostic <- look.for.trouble(param, data)
                if(!diagnostic$pass){
                    stop(paste("\n\n PARANOID MODE DETECTED AN ERROR WHILE FINDING IMPORTS:\n",
                               "at iteration ", i, ", ",
                               "movement ", j, " (", names(moves)[j], ")",
                               " with the following diagnostics:\n", diagnostic$msg, "\n\n",
                               sep=""))
                }
            }
        }

        ## store outputs if needed
        if((i %% config$sample.every) == 0 && i>1000){
            influences[counter,] <- - sapply(seq.int(data$N), function(i) ll.all(data=data, param=param, i=i))
            counter <- counter + 1L
        }

    } # end of the chain


    ## FIND OUTLIERS BASED ON INFLUENCE ##
    mean.influences <- apply(influences, 2, mean)
    mean.influence <- mean(mean.influences, na.rm=TRUE)
    threshold <- mean.influence * config$outlier.threshold
    outliers <- mean.influences > threshold


    ## RETURN ##
    ini.param$alpha[[1]][outliers] <- ini.param$current.alpha[outliers] <- NA
    ini.param$kappa[[1]][outliers] <- ini.param$current.kappa[outliers] <- NA
    config$move.alpha[outliers] <- FALSE
    config$move.kappa[outliers] <- FALSE

    return(list(config=config, param=ini.param))
} # end outbreaker.find.imports
