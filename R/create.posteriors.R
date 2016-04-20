
##
## This function creates a named list of posterior functions with enclosed likelihood and prior functions.
##

create.posteriors <- function(loglike, priors){

    ## These are all the functions generating various posterior functions;
    ## we list them by alphabetic order

    default.functions <- list(genetic = make.post.genetic,
                              reporting = make.post.reporting
                              )
    post.function.names <- names(default.functions)
    out <- lapply(default.functions, function(f) f(loglike=loglike, priors=priors))

}
