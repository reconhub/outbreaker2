
## This function is not exported and is meant for internal use in outbreaker2. It is near identical
## to the MCMC implemented by 'outbreaker.move'. The only difference is it has its own number of
## iterations ('config$n.iter.import') and sampling frequency ('config$sample.every.import'), and
## stores individual likelihoods for each case and saved iteration. The rationale is to use these
## chains to compute the 'global influences' of each case, flag outliers based on these values and
## an arbitrary threshold ('config$outlier.threshold'), and mark these cases as imported, i.e. for
## which the ancestor will be 'NA'.


outbreaker.find.imports <- function(moves, data, param, config, densities) {
    ## send back unchanged chains if config disabled the detection of imported cases ##
    if (!config$find.import) {
        return(list(config=config, param=param))
    }


    ## store initial param values ##
    ini.param <- param

    ## get number of moves ##
    J <- length(moves)

    ## create matrix of individual influences ##
    n.measures <- floor(config$n.iter.import-1000/config$sample.every.import)
    influences <- matrix(0, ncol=data$N, nrow=n.measures)
    counter <- 1L

    ## RUN MCMC ##
    for (i in seq.int(2, config$n.iter.import, 1)) {
        ## move parameters / augmented data
        for (j in seq_len(J)) {
            ## move parameters
            param <- moves[[j]](param)

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
        if ((i %% config$sample.every.import) == 0 && i>1000) {
            influences[counter,] <- - vapply(seq_len(data$N),
                                             function(i) densities$loglike$all(param=param, i=i),
                                             numeric(1))
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
}
