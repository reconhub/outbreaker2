
## This function was contributed by Rich Fitzjohn. It modifies default arguments
## using user-provided values. The argument 'strict' triggers and error
## behaviour: if strict==TRUE: all new values need to be part of the defaults.

modify_defaults <- function(defaults, x, strict = TRUE) {
    extra <- setdiff(names(x), names(defaults))
    if (strict && (length(extra) > 0L)) {
        stop("Additional invalid options: ", paste(extra, collapse=", "))
    }
    utils::modifyList(defaults, x, keep.null = TRUE) # keep.null is needed here
}






## This function pics potential infectors at random in the past, for each case.

ralpha <- function(t_inf) {
    ## choose_possible_ancestors
    canBeAnces <- outer(t_inf, t_inf, FUN = "<") # strict < is needed as we impose w(0)=0
    diag(canBeAnces) <- FALSE

    ## pick possible ancestors at random
    alpha <- apply(canBeAnces, 2,
                   function(e) ifelse(length(which(e))>0, sample(which(e),1), NA) )

    return(alpha)
}






## This function performs various checks on a MCMC, to try and pick up issues as
## soon as they appear. This includes changes to a -Inf log-likelihood or
## prior. It is meant for debugging purposes only. These checks drastically slow
## down computations.

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
        out$msg <- c(out$msg, "some delays between infections < 1 (param_current$t_inf)")
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
                                     stats::convolve(data$log_w_dens[i-1,],
                                                     rev(data$w_dens),
                                                     type="open")[seq_len(ncol(data$log_w_dens))]
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






## choose one possible ancestor for a case 'i'

.choose_possible_alpha <- function(t_inf, i) {
    return(sample(.are_possible_alpha(t_inf = t_inf, i = i), 1))
}














## ## swaps ancestries in the tree
## ## x-> i becomes i->x
## ## plus all subsequent changes

## swap_cases <- function(param, config, i) {
##     ## stop if 'i' out of range
##     if (i>length(param$alpha)) {
##         stop("trying to swap ancestry of case ",
##              i, " while there are only ",
##              length(param$alpha), " cases")
##     }

##     ## find cases for which ancestries can move
##     id_ok_to_swap <- which(can_be_swapped(param, config))

##     ## find ancestor of 'i'
##     x <- param$alpha[i]

##     ## stop if case 'i' is imported - this should not happen
##     if (is.na(x)) {
##         warning("trying to swap the ancestry of the imported case ", i)
##         return(param)
##     }

##     ## check that x can be swapped, stop if not
##     if (!(x %in% id_ok_to_swap)) {
##         return(param)
##     }

##     ## find indices to swap
##     to_be_x <- intersect(which(param$alpha==i), id_ok_to_swap)
##     to_be_i <- intersect(which(param$alpha==x), id_ok_to_swap)

##     ## swap 'i' and 'x' in ancestries
##     param$alpha[to_be_x] <- x
##     param$alpha[to_be_i] <- i

##     ## the ancestor of 'i' is now has the ancestor of 'x'
##     param$alpha[i] <- param$alpha[x]

##     ## 'i' is now the ancestor of 'x'
##     param$alpha[x] <- i

##     ## swap t_inf
##     param$t_inf[c(x,i)] <- param$t_inf[c(i,x)]

##     return(param)
## }

