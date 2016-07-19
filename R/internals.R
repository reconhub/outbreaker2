
#' Internal functions for outbreaker2
#'
#' These functions are meant for internal use in outbreaker2. They are not exported, and their API might change in the future.
#'
#' \describe{
#' \item{modify.default}{modify default arguments using user-provided values}
#' }
#'
#' @rdname internals
#'
#' @param defaults a list containing default arguments
#' @param x in \code{modify.defaults}, a list with user-provided arguments; in \code{add.to.context}, a function or an environment.
#' @param strict a logical indicating if errors shoul be returned when 'x' contains items not in 'defaults'
#'
#' @author Rich Fitzjohn, Thibaut Jombart
#'
#'
modify.defaults <- function(defaults, x, strict=TRUE) {
    extra <- setdiff(names(x), names(defaults))
    if (strict && (length(extra) > 0L)) {
        stop("Additional invalid options: ", paste(extra, collapse=", "))
    }
    modifyList(defaults, x)
}





#' @rdname internals
#' @export
#' @param t.inf a vector of infection dates
#'
ralpha <- function(t.inf) {
    ## choose.possible.ancestors
    canBeAnces <- outer(t.inf,t.inf,FUN="<") # strict < is needed as we impose w(0)=0
    diag(canBeAnces) <- FALSE

    ## pick possible ancestors at random
    alpha <- apply(canBeAnces, 2, function(e) ifelse(length(which(e))>0, sample(which(e),1), NA) )

    return(alpha)
}







##############################
## NON-DOCUMENTED FUNCTIONS ##
##############################



## checks are only sure for the 'current' state
##
look.for.trouble <- function(param.current, param.store, data) {
    ## PREPARE OUTPUT ##
    out <- list(pass=TRUE, msg=NULL)


    ## LIEKLIHOOD / POSTERIOR / PRIOR
    ## look for NAs in loglike / post / prior
    if (any(is.na(param.store$post))) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "NA detected in posterior values (param.store$post)")
    }
    if (any(is.na(param.store$like))) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "NA detected in likelihood values (param.store$like)")
    }
    if (any(is.na(param.store$prior))) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "NA detected in prior values (param.store$prior)")
    }

    ## look for NAs in loglike / post / prior
    if (!all(is.finite(param.store$post))) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "non-finite posterior values detected (param.store$post)")
    }
    if (!all(is.finite(param.store$like))) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "non-finite likelihood values detected (param.store$like)")
    }
    if (!all(is.finite(param.store$prior))) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "non-finite prior values detected (param.store$prior)")
    }


    ## CHECKS ON MU ##
    ## check that mu > 0
    if (param.current$mu<0) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "mu has a negative value:", param.current$mu)
    }

    ## check if mu is NA
    if (is.na(param.current$mu)) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "mu is NA")
    }

    ## check if mu is finite
    if (!is.finite(param.current$mu)) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "mu is not finite and equals:", param.current$mu)
    }

    ## check if mu is numeric
    if (!is.numeric(param.current$mu)) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "mu is not numeric and equals:", param.current$mu)
    }


    ## ANCESTRIES ##
    ## look for new imported cases (should not happen)
    if (!identical(is.na(param.store$alpha[[1]]), is.na(param.current$alpha))) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "imported cases have changed")
    }

    ## look for negative ancestries
    if (any(param.current$alpha<1,na.rm=TRUE)) {
       out$pass <- FALSE
       out$msg <- c(out$msg, "some ancestries point to unknown cases (param.current$alpha<1)")
    }

    ## look for ancestries greater than 'N'
    if (any(param.current$alpha>length(param.store$alpha[[1]]),na.rm=TRUE)) {
       out$pass <- FALSE
       out$msg <- c(out$msg, "some ancestries point to unknown cases (param.current$alpha>N)")
    }

    ## case infecting itself
    if (any(param.current$alpha==seq_along(param.current$alpha),na.rm=TRUE)) {
       out$pass <- FALSE
       out$msg <- c(out$msg, "auto-infections detected (param.current$alpha[i]==i)")
    }


    ## INFECTION DATES ##
    ## check NA
    if (any(is.na(param.current$t.inf))) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "NA detected in infection dates (param.current$t.inf)")
    }

    ## check finite values
    if (any(!is.finite(param.current$t.inf))) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "some infection dates are not finite (param.current$t.inf)")
    }

    ## check that values are numeric
    if (any(!is.numeric(param.current$t.inf))) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "some infection dates are not numeric (param.current$t.inf)")
    }

    ## check that delays between infections are > 0
    if (any((param.current$t.inf - param.current$t.inf[param.current$alpha]) < 1, na.rm=TRUE)) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "some delays between succesive infections are less than 1 (param.current$t.inf)")
    }

    ## check that delays to collection are > 0
    if (any((data$dates-param.current$t.inf) < 1, na.rm=TRUE)) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "some delays to collection are less than 1 (param.current$t.inf)")
    }

    ## SHAPE OUTPUT AND RETURN ##
    out$msg <- paste(out$msg, collapse="\n")
    return(out)
}





## check which ancestries can move (returns a TRUE/FALSE vector)
can.move.alpha <- function(param, config) {
    out <- !is.na(param$alpha) & # non-imported case
        (param$t.inf > min(param$t.inf)) & # not the first date
            config$move.alpha # add user-specification through move.alpha
    return(out)
}


## check which ancestries can move (returns a TRUE/FALSE vector)
can.be.swapped <- function(param, config) {
    out <- !is.na(param$alpha) & # non-imported case
            config$move.alpha # add user-specification through move.alpha
    return(out)
}


## random selection of cases for which ancestries is moved
select.alpha.to.move <- function(param, config) {
    choices <- which(can.move.alpha(param, config))
    n.to.move <- max(round(config$prop.alpha.move * length(choices)),0)
    out <- sample(choices, n.to.move, replace=FALSE)
    return(out)
}



## swaps ancestries in the tree
## x-> i becomes i->x
## plus all subsequent changes
swap.cases <- function(param, config, i) {
    ## stop if 'i' out of range
    if (i>length(param$alpha)) {
        stop("trying to swap ancestry of case ",
             i, " while there are only ",
             length(param$alpha), " cases")
    }

    ## find cases for which ancestries can move
    id.ok.to.swap <- which(can.be.swapped(param, config))

    ## find ancestor of 'i'
    x <- param$alpha[i]

    ## stop if case 'i' is imported - this should not happen
    if (is.na(x)) {
        warning("trying to swap the ancestry of the imported case ", i)
        return(param)
    }

    ## check that x can be swapped, stop if not
    if (!(x %in% id.ok.to.swap)) {
        return(param)
    }

    ## find indices to swap
    to.be.x <- intersect(which(param$alpha==i), id.ok.to.swap)
    to.be.i <- intersect(which(param$alpha==x), id.ok.to.swap)

    ## swap 'i' and 'x' in ancestries
    param$alpha[to.be.x] <- x
    param$alpha[to.be.i] <- i

    ## the ancestor of 'i' is now has the ancestor of 'x'
    param$alpha[i] <- param$alpha[x]

    ## 'i' is now the ancestor of 'x'
    param$alpha[x] <- i

    ## swap t.inf
    param$t.inf[c(x,i)] <- param$t.inf[c(i,x)]

    return(param)
}




## check that 'i' is a vector of valid case ids
## and return correct IDs
## (non-exported)
check.i <- function(data, i) {
    if (is.null(i)) seq_len(data$N) else i
    ## if (is.null(i)) return(seq_len(data$N))
    ## if (!is.numeric(i)) stop("i is not numeric")
    ## if (any(is.na(i))) stop("NA detected in case IDs")
    ## if (length(i)==0L) stop("i has length zero")
    ## if (any(i < 1)) stop("i contains invalid case indices (i<1)")
    ## if (any(i > data$N)) stop("i contains invalid case indices (i>dat$N)")
    ## return(i)
}





## find descendents of a case 'i'
find.descendents <- function(param, i) {
    ## find descendents
    which(param$alpha==i)
}





## add convolutions to data$log.w.dens
## rows = kapp avalue
## columns = time interval
add.convolutions <- function(data, config) {
    ## COMPUTE CONVOLUTIONS IF NEEDED ##
    if (config$max.kappa>1) {
        ## de-log the first line
        data$log.w.dens[1,] <- data$w.dens

        ## first compute convolutions on natural scale
        for (i in 2:config$max.kappa) {
            data$log.w.dens <- rbind(data$log.w.dens,
                                     stats::convolve(data$log.w.dens[i-1,], rev(data$w.dens), type="open")[seq_len(ncol(data$log.w.dens))]
                                     )
        }

        ## then log all densities
        data$log.w.dens <- log(data$log.w.dens)
        }

    ## name rows/columns (useful if internal debugging needed)
    rownames(data$log.w.dens) <- paste("kappa", seq_len(nrow(data$log.w.dens)), sep="=")
    colnames(data$log.w.dens) <- seq_len(ncol(data$log.w.dens))

    return(data)
}






#####################################################################################################
#####################################################################################################
## Some of the functions used in these tests have been designed in R, then translated into C++ (with
## Rcpp integration). For testing purposes, we leave the 'old' R versions here, to check the
## behaviour of new version is as it should be.
#####################################################################################################
#####################################################################################################


## This likelihood corresponds to the probability of observing a number of mutations between cases
## and their ancestors. See src/likelihoods.cpp for details of the Rcpp implmentation.

.ll.genetic <- function(data, param, i=NULL) {
    if (is.null(i)) {
        i <- seq_len(data$N)
    }

    ## discard cases with no ancestors to avoid subsetting data$D with 'NA'
    i <- i[!is.na(param$alpha[i])]

    ## likelihood is based on the number of mutations between a case and its ancestor;
    ## these are extracted from a pairwise genetic distance matrix (data$D)
    nmut <- data$D[cbind(i, param$alpha[i], deparse.level=0)]

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
    log(param$mu / (1 - param$mu)) * sum(nmut) +
        length(nmut) * log(1 - param$mu) * data$L
}






## This likelihood corresponds to the probability of observing infection dates of cases given the
## infection dates of their ancestors.

.ll.timing.infections <- function(data, param, i=NULL) {
    if (is.null(i)) {
        i <- seq_len(data$N)
    }

    ## discard cases with no ancestors to avoid subsetting data$D with 'NA'
    i <- i[!is.na(param$alpha[i])]


    ## compute delays between infection dates of cases and of their ancestors
    T <- param$t.inf[i] - param$t.inf[param$alpha[i]]

    ## avoid over-shooting: delays outside the range of columns in pre-computed log-densities
    ## (data$log.w.dens) will give a likelihood of zero
    if (any(T<1 | T>ncol(data$log.w.dens), na.rm=TRUE)) return(-Inf)

    ## output is a sum of log-densities
    sum(data$log.w.dens[cbind(param$kappa[i], T)], na.rm=TRUE)
}





## This likelihood corresponds to the probability of reporting dates of cases given their
## infection dates.

.ll.timing.sampling <- function(data, param, i=NULL) {
    if (is.null(i)) {
        i <- seq_len(data$N)
    }

    ## compute delays
    T <- data$dates[i] - param$t.inf[i]
    T <- T[!is.na(T)]

    ## avoid over-shooting
    if (any(T<1 | T>length(data$log.f.dens))) return(-Inf)

    ## output is a sum of log densities
    sum(data$log.f.dens[T], na.rm=TRUE)
}






## This likelihood corresponds to the probability of a given number of unreported cases on an ancestry.

.ll.reporting <- function(data, param, i=NULL) {
    if (is.null(i)) {
        i <- seq_len(data$N)
    }

    sum(stats::dgeom(param$kappa[i]-1,
                     prob=param$pi,
                     log=TRUE), na.rm=TRUE)
}







## This function implements movements for mu

.move.mu <- function(config, densities){
    function(param) {
        ## get new proposed values
        new.param <- param
        ##new.param$mu <- new.param$mu + rand$mu.rnorm1()
        new.param$mu <-  stats::rnorm(1, mean=new.param$mu, sd=config$sd.mu)

        ## escape if new.mu<0 or >1
        if (new.param$mu<0 || new.param$mu>1) {
            return(param)
        }

        ## compute log ratio  (assumes symmetric proposal)
        logratio <- densities$posteriors$genetic(new.param) -
            densities$posteriors$genetic(param)

        ## accept/reject
        if (logratio >= log(stats::runif(1))) {
            return(new.param)
        }
        return(param)
    }
}





## This function implements movements for t.inf

.move.t.inf <- function(config, densities) {
    prob.move <- config$prop.t.inf.move/2
    prob.proposal <- c(prob.move, 1-config$prop.t.inf.move, prob.move)

    function(param) {
        ## propose new t.inf
        new.param <- param
        new.param$t.inf <- new.param$t.inf +
            sample(-1:1, size=length(new.param$t.inf), replace=TRUE, prob=prob.proposal)

        ## compute log ratio
        logratio <- densities$loglike$timing(new.param) - densities$loglike$timing(param)

        ## accept/reject
        if (logratio >= log(stats::runif(1))) {
            return(new.param)
        } else {
            return(param)
        }
    }
}



## Movement of ancestries ('alpha') is not vectorised, movements are made one case at a time. This
## procedure is simply about picking an infector at random amongst cases preceeding the case
## considered. This movement is not symmetric, as the number of choices may change. The original
## version in 'outbreaker' used to move simultaneously 'alpha', 'kappa' and 't.inf', but current
## implementation is simpler and seems to mix at least as well. Proper movement of 'alpha' needs
## this procedure as well as a swapping procedure (swaps are not possible through move.alpha only).
##
## This is the old R function, replaced by Rcpp.move.alpha.
##
.move.alpha <- function(config, densities) {
    function(param) {
        ## create new parameters
        new.param <- param

        ## find out which ancestries to move
        alpha.can.move <- !is.na(param$alpha) & param$t.inf>min(param$t.inf)
        if (!any(alpha.can.move)) {
            warning("trying to move ancestries but none can move")
            return(param$alpha)
        }
        n.to.move <- max(round(config$prop.alpha.move * sum(alpha.can.move)),1)
        to.move <- sample(which(alpha.can.move), n.to.move, replace=FALSE)

        ## initialize new alpha
        new.param$alpha <- param$alpha

        ## move all ancestries that should be moved
        for (i in to.move) {
            ## propose new ancestor
            new.param$alpha[i] <- .choose.possible.alpha(param$t.inf, i)

            ## compute log ratio
            logratio <-  densities$loglike$all(new.param) - densities$loglike$all(param)

            ## compute correction factor
            logratio <- logratio + log(sum(.are.possible.alpha(new.param$t.inf, i))) -
                log(sum(.are.possible.alpha(param$t.inf, i)))

            ## accept/reject
            if (logratio >= log(stats::runif(1))) {
                param$alpha[i] <- new.param$alpha[i]
            } else {
                new.param$alpha[i] <- param$alpha[i]
            }
        } # end for loop

        return(param)
    }
}






## which cases are possible ancestors for a case 'i'
.are.possible.alpha <- function(t.inf, i) {
    if (length(i)>1) {
        stop("i has a length > 1")
    }
    if (any(t.inf[i]==min(t.inf))) {
        return(NA)
    }
    return(which(t.inf < t.inf[i[1]]))
}




## choose one possible ancestor for a case 'i'
.choose.possible.alpha <- function(t.inf, i) {
    return(sample(.are.possible.alpha(t.inf=t.inf, i=i), 1))
}


