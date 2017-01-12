
##
## This function creates a named list of prior functions with enclosed prior parameters
##

create_priors <- function(config) {

    ## These are all the functions generating various prior functions;
    ## we list them by alphabetic order

    default_functions <- list(mu = make_prior_mu,
                              pi = make_prior_pi
                              )
    priors_function_names <- names(default_functions)
    out <- lapply(default_functions, function(f) f(config))


    ## We add a function summing all priors

    out$all <- function(param) {
        sum(vapply(out[priors_function_names],
                   function(f) f(param), FUN.VALUE = numeric(1)), na.rm = TRUE)
    }


    return(out)
}
