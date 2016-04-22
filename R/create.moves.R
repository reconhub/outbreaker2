
##
## This function creates a named list of movement functions taking a single argument 'param'; all
## the rest (e.g. likelihood, prior, posterior functions, config, etc) is enclosed in the functions.
##

create.moves <- function(config, densities, rand){

    ## These are all the functions generating various movement functions; we list them by alphabetic
    ## order.

    ## SET DEFAULTS ##
    default.functions <- list(mu = make.move.mu,
                     t.inf = make.move.t.inf,
                     alpha = make.move.alpha,
                     swap.cases = make.move.swap.cases,
                     pi = make.move.pi,
                     kappa = make.move.kappa
                     )

    out <- lapply(default.functions, function(f) f(config, densities, rand))


    ## REMOVE FUNCTIONS IF MOVEMENTS DISABLED ##
    ## remove move$alpha if no ancestry can be moved
    if (!any(config$move.alpha)) {
        out$alpha <-  out$swap.cases <- NULL
    }

    ## remove move$t.inf if disabled
    if (!any(config$move.t.inf)) {
        out$t.inf <-  out$swap.cases <- NULL
    }

    ## remove move$mu if disabled
    if (!any(config$move.mu)) {
        out$mu <- NULL
    }

    ## remove swap if disabled, or if some alpha/t.inf cannot be moved
    if (!any(config$move.swap.cases) ||
       !any(config$move.alpha) ||
       !any(config$move.t.inf)) {
        out$swap.cases <- NULL
    }

    ## remove move$pi if disabled
    if (!any(config$move.pi)) {
        out$pi <- NULL
    }

    ## remove move$kappa if disabled
    if (!any(config$move.kappa)) {
        out$kappa <- NULL
    }


    ## return list of movement functions with enclosed stuff ##
    return(out)

}

