
##
## This function creates a named list of prior functions with enclosed prior parameters
##

create.priors <- function(config){

    ## These are all the functions generating various prior functions;
    ## we list them by alphabetic order

    default.functions <- list(mu = make.prior.mu,
                              pi = make.prior.pi
                              )
    priors.function.names <- names(default.functions)
    out <- lapply(default.functions, function(f) f(config))


    ## We add a function summing all priors

    out$all <- function(param){
        sum(vapply(out[priors.function.names], function(f) f(param), FUN.VALUE=numeric(1)), na.rm=TRUE)
    }


    return(out)
}
