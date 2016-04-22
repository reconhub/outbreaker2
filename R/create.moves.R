
##
## This function creates a named list of movement functions taking a single argument 'param'; all
## the rest (e.g. likelihood, prior, posterior functions, config, etc) is enclosed in the functions.
##

create.moves <- function(config, densities, rand){

    ## These are all the functions generating various movement functions; we list them by alphabetic
    ## order.

    ## SET DEFAULTS ##
    default.functions <- list(move.mu = make.move.mu,
                     move.t.inf = make.move.t.inf,
                     move.alpha = make.move.alpha,
                     move.swap.cases = make.move.swap.cases,
                     move.pi = make.move.pi,
                     move.kappa = make.move.kappa
                     )

    out <- lapply(default.functions, function(f) f(config, densities, rand))


    ## REMOVE FUNCTIONS IF MOVEMENTS DISABLED ##
    ## remove move.alpha if no ancestry can be moved
    if(!any(config$move.alpha)) out$move.alpha <-  out$move.swap.cases <- NULL

    ## remove move.t.inf if disabled
    if(!any(config$move.t.inf)) out$move.t.inf <-  out$move.swap.cases <- NULL

    ## remove move.mu if disabled
    if(!any(config$move.mu)) out$move.mu <- NULL

    ## remove swap if disabled, or if some alpha/t.inf cannot be moved
    if(!any(config$move.swap.cases) ||
       !any(config$move.alpha) ||
       !any(config$move.t.inf)) out$move.swap.cases <- NULL

    ## remove move.pi if disabled
    if(!any(config$move.pi)) out$move.pi <- NULL

    ## remove move.kappa if disabled
    if(!any(config$move.kappa)) out$move.kappa <- NULL


    ## return list of movement functions with enclosed stuff ##
    return(out)

}

