
##
## This function creates a named list of posterior functions with enclosed likelihood and prior functions.
##

create.posteriors <- function(loglike, priors){

    ## These are all the functions generating various posterior functions; we list them by
    ## alphabetic order; note that unlike likelihood and prior functions, the '$all' here is not the
    ## sum of the components of the list, but directly uses the 'all' functions in 'loglike' and
    ## 'priors'.


    default.functions <- list(genetic = make.post.genetic,
                              reporting = make.post.reporting,
                              all = make.post.all)

    out <- lapply(default.functions, function(f) f(loglike=loglike, priors=priors))

}
