
##
## This function creates a named list of movement functions taking a single argument 'param'; all
## the rest (e.g. likelihood, prior, posterior functions, config, etc) is enclosed in the functions.
##

create.moves <- function(config, densities, rand){

    ## These are all the functions generating various movement functions; we list them by alphabetic
    ## order.

    ## SET DEFAULTS ##
    default.functions <- list(move.mu = move.mu,
                     move.t.inf = move.t.inf,
                     move.alpha = move.alpha,
                     move.swap.cases = move.swap.cases,
                     move.pi = move.pi,
                     move.kappa = move.kappa
                     )

    out <- lapply(default.functions, function(f) f(config, densities, rand))


    ## REMOVE FUNCTIONS IF MOVEMENTS DISABLED ##
    ## remove move.alpha if no ancestry can be moved
    if(!any(config$move.alpha)) moves$move.alpha <-  moves$move.swap.cases <- NULL

    ## remove move.t.inf if disabled
    if(!any(config$move.t.inf)) moves$move.t.inf <-  moves$move.swap.cases <- NULL

    ## remove move.mu if disabled
    if(!any(config$move.mu)) moves$move.mu <- NULL

    ## remove swap if disabled, or if some alpha/t.inf cannot be moved
    if(!any(config$move.swap.cases) ||
       !any(config$move.alpha) ||
       !any(config$move.t.inf)) moves$move.swap.cases <- NULL

    ## remove move.pi if disabled
    if(!any(config$move.pi)) moves$move.pi <- NULL

    ## remove move.kappa if disabled
    if(!any(config$move.kappa)) moves$move.kappa <- NULL


    ## ## CHECK FUNCTIONS ##
    ## check.function.args <- function(f){
    ##     if(is.null(f)) return(TRUE)
    ##     args <- names(formals(f))
    ##     if(!identical(sort(args), c("config","data","param", "rand"))) {
    ##         return(FALSE)
    ##     } else {
    ##         return(TRUE)
    ##     }
    ## }

    ## args.checks <- sapply(unlist(moves), check.function.args)
    ## if(!all(args.checks)){
    ##     culprits <- names(args.checks)[!args.checks]
    ##     culprits <- paste(culprits, collapse=",")
    ##     stop("problems in movements of: ", culprits, "\narguments shoud be: 'data', 'config', 'param', 'rand'")
    ## }

    ## RETURN ##
    return(moves)

} # end outbreaker.create.moves

