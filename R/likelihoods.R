

## In all of these functions:
## --------------------------
##
## we return a function in which 'data' is enclosed this is substantially faster than having to
## pass the extra 'data' argument all the time, which does not change anyway; only things which
## effectively change are the parameters and augmented data, all stored in 'param', and the
## indices of cases for which we compute the likelihood 'i'.

## Arguments are:
##
## 'data' a list of named items containing input data as returned by \code{\link{outbreaker.data}}
## 'param' a list containing parameters as returned by \code{outbreaker.create.mcmc}
## 'i' an optional vector of integers, indicating subset of cases included in the likelihood computation; if NULL (default), all cases are used.




## This likelihood corresponds to the probability of observing infection dates of cases given the
## infection dates of their ancestors.

make.ll.timing.infections <- function(data) {
    ## i will be the index of cases to be used, but it is useful to define it by default as all cases
    cases <- seq_len(data$N)

    if (data$N>1) {
        function(param, i=cases) {

            ## compute delays between infection dates of cases and of their ancestors
            T <- param$current.t.inf[i] - param$current.t.inf[param$current.alpha[i]]

            ## avoid over-shooting: delays outside the range of columns in pre-computed log-densities
            ## (data$log.w.dens) will give a likelihood of zero
            if (any(T<1 | T>ncol(data$log.w.dens), na.rm=TRUE)) return(-Inf)

            ## return sum of log-densities
            sum(data$log.w.dens[cbind(param$current.kappa[i], T)], na.rm=TRUE)
        }
    } else {
        function(...) 0
    }
}




## This likelihood corresponds to the probability of reporting dates of cases given their
## infection dates.

make.ll.timing.sampling <- function(data) {
    ## i will be the index of cases to be used, but it is useful to define it by default as all cases
    cases <- seq_len(data$N)

    if (data$N>1) {
        function(param, i=cases) {

            ## compute delays
            T <- data$dates[i] - param$current.t.inf[i]
            T <- T[!is.na(T)]

            ## avoid over-shooting
            if (any(T<1 | T>length(data$log.f.dens))) return(-Inf)

            ## return
            sum(data$log.f.dens[T], na.rm=TRUE)
        }
    } else {
        function(...) 0
    }
}





## This likelihood corresponds to the probability of observing a number of mutations between cases
## and their ancestors.

make.ll.genetic <- function(data) {
    ## i will be the index of cases to be used, but it is useful to define it by default as all cases
    cases <- seq_len(data$N)

    if (nrow(data$D)>1) {

        function(param, i=cases) {

            ## discard cases with no ancestors to avoid subsetting data$D with 'NA'
            i <- i[!is.na(param$current.alpha[i])]

            ## likelihood is based on the number of mutations between a case and its ancestor;
            ## these are extracted from a pairwise genetic distance matrix (data$D)
            nmut <- data$D[cbind(i, param$current.alpha[i], deparse.level=0)]

            ## the log-likelihood is computed as: sum(mu^nmut + (1-mu)^(L-nmut))
            ## with:
            ## 'mu' is the mutation probability
            ## 'L' the number of sites in the alignment
            ## 'nmut' the number of mutations between an ancestor and its descendent
            ##
            ## for computer efficiency, we re-factorise it as:
            ##  log(mu / (1 - mu)) * sum(nmut) + length(nmut) * log(1 - mu) * L
            ## which limits to 2 operations rather than 2*n
            ## (tip from Rich Fitzjohn)
            log(param$current.mu / (1 - param$current.mu)) * sum(nmut) + length(nmut) * log(1 - param$current.mu) * data$L
        }
    } else {
        function(...) 0
    }
}





## This likelihood corresponds to the probability of a given number of unreported cases on an ancestry.

make.ll.reporting <- function(data) {
    ## i will be the index of cases to be used, but it is useful to define it by default as all cases
    cases <- seq_len(data$N)

    ## the likelihood is given by a geometric distribution with probability 'pi' to report a case
    ## 'kappa' is the number of generation between two successive cases
    ## 'kappa-1' is the number of unreported cases
    if (data$N>1) {

        function(param, i=cases) {
            sum(stats::dgeom(param$current.kappa[i]-1,
                             prob=param$current.pi,
                             log=TRUE), na.rm=TRUE)
        }
    } else {
        function(...) 0
    }
}
