
#' Internal functions for outbreaker2
#'
#' These functions are meant for internal use in outbreaker2. They are not
#' exported, and their API might change in the future.
#'
#' \describe{
#' \item{modify_default}{modify default arguments using user-provided values}
#' }
#'
#' @rdname internals
#'
#' @param defaults a list containing default arguments
#' 
#' @param x in \code{modify_defaults}, a list with user-provided arguments; in
#'     \code{add_to_context}, a function or an environment.
#'
#' @param strict a logical indicating if errors shoul be returned when 'x'
#'     contains items not in 'defaults'
#'
#' @author Rich Fitzjohn, Thibaut Jombart
#'
#'
modify_defaults <- function(defaults, x, strict = TRUE) {
    extra <- setdiff(names(x), names(defaults))
    if (strict && (length(extra) > 0L)) {
        stop("Additional invalid options: ", paste(extra, collapse=", "))
    }
    modifyList(defaults, x)
}





#' @rdname internals
#' @export
#' @param t_inf a vector of infection dates
#'
ralpha <- function(t_inf) {
    ## choose_possible_ancestors
    canBeAnces <- outer(t_inf, t_inf, FUN = "<") # strict < is needed as we impose w(0)=0
    diag(canBeAnces) <- FALSE

    ## pick possible ancestors at random
    alpha <- apply(canBeAnces, 2,
                   function(e) ifelse(length(which(e))>0, sample(which(e),1), NA) )

    return(alpha)
}







##############################
## NON-DOCUMENTED FUNCTIONS ##
##############################



## checks are only sure for the 'current' state
##
look_for_trouble <- function(param_current, param_store, data) {
    ## PREPARE OUTPUT ##
    out <- list(pass = TRUE, msg = NULL)


    ## LIEKLIHOOD / POSTERIOR / PRIOR
    ## look for NAs in loglike / post / prior
    if (any(is.na(param_store$post))) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "NA detected in posterior values (param_store$post)")
    }
    if (any(is.na(param_store$like))) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "NA detected in likelihood values (param_store$like)")
    }
    if (any(is.na(param_store$prior))) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "NA detected in prior values (param_store$prior)")
    }

    ## look for NAs in loglike / post / prior
    if (!all(is.finite(param_store$post))) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "non-finite posterior values detected (param_store$post)")
    }
    if (!all(is.finite(param_store$like))) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "non-finite likelihood values detected (param_store$like)")
    }
    if (!all(is.finite(param_store$prior))) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "non-finite prior values detected (param_store$prior)")
    }


    ## CHECKS ON MU ##
    ## check that mu > 0
    if (param_current$mu<0) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "mu has a negative value:", param_current$mu)
    }

    ## check if mu is NA
    if (is.na(param_current$mu)) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "mu is NA")
    }

    ## check if mu is finite
    if (!is.finite(param_current$mu)) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "mu is not finite and equals:", param_current$mu)
    }

    ## check if mu is numeric
    if (!is.numeric(param_current$mu)) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "mu is not numeric and equals:", param_current$mu)
    }


    ## ANCESTRIES ##
    ## look for new imported cases (should not happen)
    if (!identical(is.na(param_store$alpha[[1]]), is.na(param_current$alpha))) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "imported cases have changed")
    }

    ## look for negative ancestries
    if (any(param_current$alpha<1,na.rm = TRUE)) {
       out$pass <- FALSE
       out$msg <- c(out$msg, "some ancestries point to unknown cases (param_current$alpha<1)")
    }

    ## look for ancestries greater than 'N'
    if (any(param_current$alpha>length(param_store$alpha[[1]]),na.rm = TRUE)) {
       out$pass <- FALSE
       out$msg <- c(out$msg, "some ancestries point to unknown cases (param_current$alpha>N)")
    }

    ## case infecting itself
    if (any(param_current$alpha==seq_along(param_current$alpha),na.rm = TRUE)) {
       out$pass <- FALSE
       out$msg <- c(out$msg, "auto-infections detected (param_current$alpha[i]==i)")
    }


    ## INFECTION DATES ##
    ## check NA
    if (any(is.na(param_current$t_inf))) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "NA detected in infection dates (param_current$t_inf)")
    }

    ## check finite values
    if (any(!is.finite(param_current$t_inf))) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "some infection dates are not finite (param_current$t_inf)")
    }

    ## check that values are numeric
    if (any(!is.numeric(param_current$t_inf))) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "some infection dates are not numeric (param_current$t_inf)")
    }

    ## check that delays between infections are > 0
    if (any((param_current$t_inf - param_current$t_inf[param_current$alpha]) < 1, na.rm = TRUE)) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "some delays between succesive infections are less than 1 (param_current$t_inf)")
    }

    ## check that delays to collection are > 0
    if (any((data$dates-param_current$t_inf) < 1, na.rm = TRUE)) {
        out$pass <- FALSE
        out$msg <- c(out$msg, "some delays to collection are less than 1 (param_current$t_inf)")
    }

    ## SHAPE OUTPUT AND RETURN ##
    out$msg <- paste(out$msg, collapse="\n")
    return(out)
}





## check which ancestries can move (returns a TRUE/FALSE vector)
can_move_alpha <- function(param, config) {
    out <- !is.na(param$alpha) & # non-imported case
        (param$t_inf > min(param$t_inf)) & # not the first date
            config$move_alpha # add user-specification through move_alpha
    return(out)
}


## check which ancestries can move (returns a TRUE/FALSE vector)
can_be_swapped <- function(param, config) {
    out <- !is.na(param$alpha) & # non-imported case
            config$move_alpha # add user-specification through move_alpha
    return(out)
}


## random selection of cases for which ancestries is moved
select_alpha_to_move <- function(param, config) {
    choices <- which(can_move_alpha(param, config))
    n_to_move <- max(round(config$prop_alpha_move * length(choices)),0)
    out <- sample(choices, n_to_move, replace = FALSE)
    return(out)
}



## swaps ancestries in the tree
## x-> i becomes i->x
## plus all subsequent changes
swap_cases <- function(param, config, i) {
    ## stop if 'i' out of range
    if (i>length(param$alpha)) {
        stop("trying to swap ancestry of case ",
             i, " while there are only ",
             length(param$alpha), " cases")
    }

    ## find cases for which ancestries can move
    id_ok_to_swap <- which(can_be_swapped(param, config))

    ## find ancestor of 'i'
    x <- param$alpha[i]

    ## stop if case 'i' is imported - this should not happen
    if (is.na(x)) {
        warning("trying to swap the ancestry of the imported case ", i)
        return(param)
    }

    ## check that x can be swapped, stop if not
    if (!(x %in% id_ok_to_swap)) {
        return(param)
    }

    ## find indices to swap
    to_be_x <- intersect(which(param$alpha==i), id_ok_to_swap)
    to_be_i <- intersect(which(param$alpha==x), id_ok_to_swap)

    ## swap 'i' and 'x' in ancestries
    param$alpha[to_be_x] <- x
    param$alpha[to_be_i] <- i

    ## the ancestor of 'i' is now has the ancestor of 'x'
    param$alpha[i] <- param$alpha[x]

    ## 'i' is now the ancestor of 'x'
    param$alpha[x] <- i

    ## swap t_inf
    param$t_inf[c(x,i)] <- param$t_inf[c(i,x)]

    return(param)
}




## check that 'i' is a vector of valid case ids
## and return correct IDs
## (non-exported)
check_i <- function(data, i) {
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
find_descendents <- function(param, i) {
    ## find descendents
    which(param$alpha==i)
}





## add convolutions to data$log_w_dens
## rows = kapp avalue
## columns = time interval
add_convolutions <- function(data, config) {
    ## COMPUTE CONVOLUTIONS IF NEEDED ##
    if (config$max_kappa>1) {
        ## de-log the first line
        data$log_w_dens[1,] <- data$w_dens

        ## first compute convolutions on natural scale
        for (i in 2:config$max_kappa) {
            data$log_w_dens <- rbind(data$log_w_dens,
                                     stats::convolve(data$log_w_dens[i-1,], rev(data$w_dens), type="open")[seq_len(ncol(data$log_w_dens))]
                                     )
        }

        ## then log all densities
        data$log_w_dens <- log(data$log_w_dens)
        }

    ## name rows/columns (useful if internal debugging needed)
    rownames(data$log_w_dens) <- paste("kappa", seq_len(nrow(data$log_w_dens)), sep="=")
    colnames(data$log_w_dens) <- seq_len(ncol(data$log_w_dens))

    return(data)
}






################################################################################
## Some of the functions used in these tests have been designed in R, then
## translated into C++ (with Rcpp integration). For testing purposes, we leave
## the 'old' R versions here, to check the behaviour of new version is as it
## should be.
################################################################################


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






## This likelihood corresponds to the probability of a given number of unreported cases on an ancestry.

.ll_reporting <- function(data, param, i = NULL) {
    if (is.null(i)) {
        i <- seq_len(data$N)
    }

    sum(stats::dgeom(param$kappa[i]-1,
                     prob = param$pi,
                     log = TRUE), na.rm = TRUE)
}







## This function implements movements for mu

.move_mu <- function(config, densities){
    function(param) {
        ## get new proposed values
        new_param <- param
        ##new_param$mu <- new_param$mu + rand$mu_rnorm1()
        new_param$mu <-  stats::rnorm(1, mean = new_param$mu, sd = config$sd_mu)

        ## escape if new_mu<0 or >1
        if (new_param$mu<0 || new_param$mu>1) {
            return(param)
        }

        ## compute log ratio  (assumes symmetric proposal)
        logratio <- densities$posteriors$genetic(new_param) -
            densities$posteriors$genetic(param)

        ## accept/reject
        if (logratio >= log(stats::runif(1))) {
            return(new_param)
        }
        return(param)
    }
}





## This function implements movements for pi.

.move_pi <- function(config, densities){
    function(param) {
        ## get new proposed values
        new_param <- param
        ## new_param$pi <- new_param$pi + rand$pi_rnorm1()
        new_param$pi <- stats::rnorm(1, mean = new_param$pi, sd = config$sd_pi)

        ## escape if new_pi<0 or >1
        if (new_param$pi<0 || new_param$pi>1) {
            return(param)
        }

        ## compute log ratio  (assumes symmetric proposal)
        logratio <- densities$posteriors$reporting(new_param) -
            densities$posteriors$reporting(param)

        ## accept/reject
        if (logratio >= log(stats::runif(1))) {
            return(new_param)
        }
        return(param)
    }
}





## This function implements movements for t_inf

.move_t_inf <- function(config, densities) {
    prob_move <- config$prop_t_inf_move/2
    prob_proposal <- c(prob_move, 1-config$prop_t_inf_move, prob_move)

    function(param) {
        ## propose new t_inf
        new_param <- param
        new_param$t_inf <- new_param$t_inf +
            sample(-1:1, size = length(new_param$t_inf), replace = TRUE, prob = prob_proposal)

        ## compute log ratio
        logratio <- densities$loglike$timing(new_param) - densities$loglike$timing(param)

        ## accept/reject
        if (logratio >= log(stats::runif(1))) {
            return(new_param)
        } else {
            return(param)
        }
    }
}



## Movement of ancestries ('alpha') is not vectorised, movements are made one
## case at a time. This procedure is simply about picking an infector at random
## amongst cases preceeding the case considered. This movement is not symmetric,
## as the number of choices may change. The original version in 'outbreaker'
## used to move simultaneously 'alpha', 'kappa' and 't_inf', but current
## implementation is simpler and seems to mix at least as well. Proper movement
## of 'alpha' needs this procedure as well as a swapping procedure (swaps are
## not possible through move_alpha only).
##
## This is the old R function, replaced by Rcpp_move_alpha.
##
.move_alpha <- function(config, densities) {
    function(param) {
        ## create new parameters
        new_param <- param

        ## find out which ancestries to move
        alpha_can_move <- !is.na(param$alpha) & param$t_inf>min(param$t_inf)
        if (!any(alpha_can_move)) {
            warning("trying to move ancestries but none can move")
            return(param$alpha)
        }
        n_to_move <- max(round(config$prop_alpha_move * sum(alpha_can_move)), 1)
        to_move <- sample(which(alpha_can_move), n_to_move, replace = FALSE)

        ## initialize new alpha
        new_param$alpha <- param$alpha

        ## move all ancestries that should be moved
        for (i in to_move) {
            ## propose new ancestor
            new_param$alpha[i] <- .choose_possible_alpha(param$t_inf, i)

            ## compute log ratio
            logratio <-  densities$loglike$all(new_param) - densities$loglike$all(param)

            ## compute correction factor
            logratio <- logratio + log(sum(.are_possible_alpha(new_param$t_inf, i))) -
                log(sum(.are_possible_alpha(param$t_inf, i)))

            ## accept/reject
            if (logratio >= log(stats::runif(1))) {
                param$alpha[i] <- new_param$alpha[i]
            } else {
                new_param$alpha[i] <- param$alpha[i]
            }
        } # end for loop

        return(param)
    }
}





## This function implements moves for kappa.

.move_kappa <- function(config, densities) {
    function(param) {
        ## browser()
        ## determine which cases to move
        kappa_can_move <- !is.na(param$kappa)
        n_to_move <- max(round(.2 * sum(kappa_can_move)), 1)
        to_move <- sample(which(kappa_can_move), n_to_move, replace = FALSE)

        ## initialize new kappa
        new_param <- param

        ## move all ancestries that should be moved
        for (i in to_move) {
            ## propose new kappa
            new_param$kappa[i] <- new_param$kappa[i] + sample(c(-1,1), size = 1)

            ## reject move automatically if new kappa < 1 or greater than allowed max
            if (new_param$kappa[i] < 1 ||
                new_param$kappa[i] > config$max_kappa) {
                new_param$kappa[i] <- param$kappa[i]
            } else {
                ## compute log ratio
                logratio <- densities$loglike$timing_infections(new_param, i) +
                    densities$loglike$genetic(new_param, i) +
                    densities$loglike$reporting(new_param, i) -
                    densities$loglike$timing_infections(param, i) -
                    densities$loglike$genetic(param, i) -
                    densities$loglike$reporting(param, i)

                ## accept/reject
                if (logratio >= log(stats::runif(1))) {
                    param$kappa[i] <- new_param$kappa[i]
                } else {
                    new_param$kappa[i] <- param$kappa[i]
                }
            }
        } # end for loop


        ## output is a list of potentially modified parameters
        return(param)
    }
}






## which cases are possible ancestors for a case 'i'
.are_possible_alpha <- function(t_inf, i) {
    if (length(i)>1) {
        stop("i has a length > 1")
    }
    if (any(t_inf[i]==min(t_inf))) {
        return(NA)
    }
    return(which(t_inf < t_inf[i[1]]))
}




## This is the complementary procedure to the above one (move_alpha). This type
## of move swaps a case 'a' with its ancestor, e_g.

## x -> a -> b  becomes a -> x -> b

## Obviously cases are moved one at a time. We need to used local likelihood
## changes for this move to scale well with outbreak size. The complicated bit
## is that the move impacts all descendents from 'a' as well as 'x'.

.move_swap_cases <- function(config, densities) {
    function(param) {
    ## find ancestries which can move
    to_move <- select_alpha_to_move(param, config)

    ## leave if nothing moves
    if (length(to_move)<1) {
        return(param)
    }

    ## move all ancestries that should be moved
    for (i in to_move) {
        ## swap ancestries
        new_param <- swap_cases(param, config, i)

        ## compute log ratio using local changes only; these include:

        ## descendents of to_move
        ## descendents of alpha[to_move]
        ## alpha[to_move]

        affected_cases <- c(find_descendents(param, i = i),
                            find_descendents(param, i = param$alpha[i]),
                            param$alpha[i])
        logratio <- densities$loglike$all(new_param, i = affected_cases) -
            densities$loglike$all(param, i = affected_cases)

        ## accept/reject
        if (logratio >= log(stats::runif(1))) {
            param <- new_param
        }
    } # end for loop

    return(param)
    }
}






## choose one possible ancestor for a case 'i'
.choose_possible_alpha <- function(t_inf, i) {
    return(sample(.are_possible_alpha(t_inf = t_inf, i = i), 1))
}
