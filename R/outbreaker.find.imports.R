#' Detection of imported cases for outbreaker2
#'
#' This function moves all parameters for outbreaker2.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @param moves a set of movement functions stored in a named list as returned by \code{\link{outbreaker.create.moves}}
#' @param data a list of data items as returned by \code{outbreaker.data}
#' @param param a list of parameters as returned by \code{outbreaker.create.mcmc}
#' @param config a list of settings as returned by \code{outbreaker.config}
#'
#' @return a potentially modified list of parameters as returned by \code{outbreaker.create.mcmc}
#'
outbreaker.find.imports <- function(moves, data, param, config, densities){
    ## send back unchanged chains if config disabled the detection of imported cases ##
    if (!config$find.import) return(list(config=config, param=param))


    ## store initial param values ##
    ini.param <- param

    ## get number of moves ##
    J <- length(moves)

    ## create matrix of individual influences ##
    n.measures <- floor(config$n.iter.import-1000/config$sample.every.import)
    influences <- matrix(0, ncol=data$N, nrow=n.measures)
    counter <- 1L

    ## RUN MCMC ##
    for(i in seq.int(2, config$n.iter.import, 1)){
        ## move parameters / augmented data
        for(j in seq.int(J)){
            ## safemode
            if (config$paranoid){
                old.param <- param
            }

            ## move parameters
            param <- moves[[j]](param)

            ## safemode
            if (config$paranoid){
                diagnostic <- look.for.trouble(param, data)
                if (!diagnostic$pass){
                    stop(paste("\n\n PARANOID MODE DETECTED AN ERROR WHILE FINDING IMPORTS:\n",
                               "at iteration ", i, ", ",
                               "movement ", j, " (", names(moves)[j], ")",
                               " with the following diagnostics:\n", diagnostic$msg, "\n\n",
                               sep=""))
                }
            }
        }

        ## store outputs if needed
        if ((i %% config$sample.every.import) == 0 && i>1000){
            influences[counter,] <- - sapply(seq.int(data$N),
                                             function(i) densities$loglike$all(param=param, i=i))
            counter <- counter + 1L
        }

    } # end of the chain


    ## FIND OUTLIERS BASED ON INFLUENCE ##
    mean.influences <- apply(influences, 2, mean)
    mean.influence <- mean(mean.influences, na.rm=TRUE)
    threshold <- mean.influence * config$outlier.threshold
    outliers <- mean.influences > threshold


    ## All outliers are considered as introductions, so that ancestries (alpha) are set to 'NA' and
    ## the number of generations between cases and their ancestor (kappa) is set to NA; the
    ## movements of alpha and kappa for these cases is also disabled; because the config has been
    ## altered in these cases, we systematically return the config as well as the initial
    ## parameters.

    ini.param$alpha[[1]][outliers] <- ini.param$current.alpha[outliers] <- NA
    ini.param$kappa[[1]][outliers] <- ini.param$current.kappa[outliers] <- NA
    ini.param$influences <- mean.influence
    ini.param$threshold <- threshold
    config$move.alpha[outliers] <- FALSE
    config$move.kappa[outliers] <- FALSE

    return(list(config=config, param=ini.param))
} # end outbreaker.find.imports
