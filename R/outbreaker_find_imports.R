
## This function is not exported and is meant for internal use in
## outbreaker2. It is near identical to the MCMC implemented by
## 'outbreaker.move'. The only difference is it has its own number of iterations
## ('config$n.iter.import') and sampling frequency
## ('config$sample.every.import'), and stores individual likelihoods for each
## case and saved iteration. The rationale is to use these chains to compute the
## 'global influences' of each case, flag outliers based on these values and an
## arbitrary threshold ('config$outlier.threshold'), and mark these cases as
## imported, i.e. for which the ancestor will be 'NA'.


outbreaker_find_imports <- function(moves, data, param_current,
                                    param_store, config,
                                    likelihoods) {
  ## send back unchanged chains if config disabled the detection of imported cases ##
  if (!config$find_import) {
    return(list(config = config,
                param_current = param_current,
                param_store = param_store))
  }

  ## store initial param values ##
  ini_param <- list(current = param_current, store = param_store)

  ## get number of moves ##
  J <- length(moves)

  ## create matrix of individual influences ##
  n_measures <- floor((config$n_iter_import - 1000) / config$sample_every_import)
  influences <- matrix(0, ncol = data$N, nrow = n_measures)
  counter <- 1L

  ## remove the contact likelihood from outlier detection
  tmp_likelihoods <- likelihoods
  tmp_likelihoods$contact <- function(data, param) return(0)

  ## remove the hospital transfer likelihood from outlier detection
  tmp_likelihoods$patient_transfer <- function(data, param) return(0)

  ## RUN MCMC ##
  for (i in seq.int(2, config$n_iter_import, 1)) {
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
    if ((i %% config$sample_every_import) == 0 && i>1000) {
      influences[counter,] <- - vapply(seq_len(data$N),
                                       function(i) cpp_ll_all(data, param_current, i, tmp_likelihoods),
                                       numeric(1))
      counter <- counter + 1L
    }
  } # end of the chain

  ## FIND OUTLIERS BASED ON INFLUENCE ##
  mean_influences <- colMeans(influences)
  mean_influence <- mean(mean_influences, na.rm = TRUE)
  threshold <- mean_influence * config$outlier_threshold
  outliers <- mean_influences > threshold


  ## All outliers are considered as introductions, so that ancestries (alpha) are set to 'NA' and
  ## the number of generations between cases and their ancestor (kappa) is set to NA; the
  ## movements of alpha and kappa for these cases is also disabled; because the config has been
  ## altered in these cases, we systematically return the config as well as the initial
  ## parameters.

  ini_param$store$alpha[[1]][outliers] <- ini_param$current$alpha[outliers] <- NA
  ini_param$store$kappa[[1]][outliers] <- ini_param$current$kappa[outliers] <- NA
  ini_param$store$influences <- mean_influence
  ini_param$store$threshold <- threshold
  config$move_alpha[outliers] <- FALSE
  config$move_kappa[outliers] <- FALSE

  return(list(config = config,
              param_current = ini_param$current,
              param_store = ini_param$store))
}
