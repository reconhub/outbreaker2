
## ON THE PURPOSE OF THESE FUNCTIONS

## These functions are no longer used in outbreaker2, but were part of the
## original implementation, and are still used in testing procedures to ensure
## that Rcpp versions give identical results.






## This likelihood corresponds to the probability of observing a number of
## mutations between cases and their ancestors. See src/likelihoods_cpp for
## details of the Rcpp implmentation.

.ll_genetic <- function(data, param, i = NULL) {
    if (is.null(i)) {
        i <- seq_len(data$N)
    }

    ## discard cases with no ancestors to avoid subsetting data$D with 'NA'
    i <- i[!is.na(param$alpha[i])]

    ## likelihood is based on the number of mutations between a case and its
    ## ancestor; these are extracted from a pairwise genetic distance matrix
    ## (data$D)
    
    ## the log-likelihood is computed as: sum(mu^n_mut + (1-mu)^(L-n_mut))
    ## with:
    ## 'mu' is the mutation probability
    ## 'L' the number of sites in the alignment
    ## 'n_mut' the number of mutations between an ancestor and its descendent
    ##
    ## for computer efficiency, we re-factorise it as:
    ##  log(mu) * sum(n_mut) + log(1 - mu) * (L-n_mut + (kappa-1)*L)
    ##

    n_mut <- data$D[cbind(i, param$alpha[i], deparse.level = 0)]
    
    n_non_mut <- (data$L - n_mut) + (param$kappa[i]-1) * data$L

    out <- log(param$mu) * sum(n_mut, na.rm = TRUE) + # mutated sites
        log(1 - param$mu) * sum(n_non_mut, na.rm = TRUE) # non mutated sites
    return(out)
}






## This likelihood corresponds to the probability of observing infection dates of cases given the
## infection dates of their ancestors.

.ll_timing_infections <- function(data, param, i = NULL) {
    if (is.null(i)) {
        i <- seq_len(data$N)
    }

    ## discard cases with no ancestors to avoid subsetting data$D with 'NA'
    i <- i[!is.na(param$alpha[i])]


    ## compute delays between infection dates of cases and of their ancestors
    T <- param$t_inf[i] - param$t_inf[param$alpha[i]]

    ## avoid over-shooting: delays outside the range of columns in pre-computed log-densities
    ## (data$log_w_dens) will give a likelihood of zero
    if (any(T<1 | T>ncol(data$log_w_dens), na.rm = TRUE)) return(-Inf)

    ## output is a sum of log-densities
    sum(data$log_w_dens[cbind(param$kappa[i], T)], na.rm = TRUE)
}





## This likelihood corresponds to the probability of reporting dates of cases given their
## infection dates.

.ll_timing_sampling <- function(data, param, i = NULL) {
    if (is.null(i)) {
        i <- seq_len(data$N)
    }

    ## compute delays
    T <- data$dates[i] - param$t_inf[i]
    T <- T[!is.na(T)]

    ## avoid over-shooting
    if (any(T<1 | T>length(data$log_f_dens))) return(-Inf)

    ## output is a sum of log densities
    sum(data$log_f_dens[T], na.rm = TRUE)
}






## This likelihood corresponds to the probability of a given number of
## unreported cases on an ancestry.

.ll_reporting <- function(data, param, i = NULL) {
    if (is.null(i)) {
        i <- seq_len(data$N)
    }

    sum(stats::dgeom(param$kappa[i]-1,
                     prob = param$pi,
                     log = TRUE), na.rm = TRUE)
}


