

## This function creates a named list of prior functions with enclosed prior
## parameters. It accepts an argument ... which can be used to pass custom prior
## functions. If not provided, then default priors are built from using
## 'config'. If provided, function must have a single argument 'param', which is
## a list containing all the parameters of the models. Parameters of the priors
## need to be either hard-coded or enclosed.


create_priors <- function(config, ...) {

    dots <- list(...)
    out <- list(...)

    ## checks: make sure all the ... are functions and they have a single
    ## argument '...'

    if (!is.null(dots)) {
        
        ## check all objects are functions
        all.ok <- all(vapply(dots, is.function, logical(1)))
        if (!all.ok) {
            msg <- "Some of the provided custom priors are not functions."
            stop("msg")
        }

        ## check they all have a single parameter 'param'
        all.ok <- all(vapply(dots, formalArgs, character(1)) == "param")
        if (!all.ok) {
            msg <- paste("All custom priors should have a single argument",
                         "'param'.")
            stop("msg")
        }
    }
    

    ## prior for mutation rate
    if (is.null(dots$mu)) {
        out$mu <- bind_prior_mu(config)
    }


    ## prior for reporting rate
    if (is.null(dots$pi)) {
        out$pi <- bind_prior_pi(config)
    }

    priors_function_names <- names(out)


    ## Function summing all priors

    out$all <- function(param) {
        sum(vapply(out[priors_function_names],
                   function(f) f(param), FUN.VALUE = numeric(1)), na.rm = TRUE)
    }


    return(out)
}
