

## This function creates a named list of prior functions with enclosed prior
## parameters. It accepts an argument ... which can be used to pass custom prior
## functions. If not provided, then default priors are built from using
## 'config'. If provided, function must have a single argument 'param', which is
## a list containing all the parameters of the models. Parameters of the priors
## need to be either hard-coded or enclosed.


create_priors <- function(config = NULL, ...) {

    ## This function returns a list of functions with the class
    ## 'outbreaker_priors'. It optionally generates priors and tests some basic
    ## properties of the prior functions. There are three ways a user can
    ## specify priors:

    ## 1) Defaults: this is what happens when the 'config' has default values of
    ## prior parameters.

    ## 2) Customized parameters: in this case, the prior functions are the
    ## default ones from the package, but will use custom parameters, passed
    ## through config$prior_[parameter name]

    ## 3) Customized functions: in this case, prior functions themselves are
    ## specified by the user, through the '...' argument. The requirements is
    ## that such functions must have either hard-coded parameters or enclosed
    ## values. They will take a single argument 'param', which is a list
    ## containing all model parameters with the class 'outbreaker_param'.


    ## Get user-specified prior functions
    
    priors <- list(...)
    if (length(priors) == 1L && is.list(priors[[1]])) {
        priors <- priors[[1]]
    }
    

    ## Handle config; if config is NULL this will create a default config;
    ## otherwise it processes it and goes through a bunch of checks

    config <- create_config(config)


    ## Use user-provided priors where provided, default otherwise. The default
    ## for a prior is NULL, in which case the movement functions in C++ will use
    ## C++ versions.
    
    defaults <- list(mu = NULL, # mutation rate
                     pi = NULL # reporting probability
                     )
    
    priors <- modify_defaults(defaults, priors, FALSE)
    priors_names <- setdiff(names(priors), "all")

    ## ## prior for mutation rate

    ## if (is.null(priors$mu)) {
    ##     priors$mu <- bind_prior_mu(config)
    ## }

    
    ## ## prior for reporting rate
    
    ## if (is.null(priors$pi)) {
    ##     priors$pi <- bind_prior_pi(config)
    ## }

    
    ## Run checks; these include making sure that:

    ## - all priors are there
    ## - all priors are functions
    ## - all priors have a single argument 'param'


    ## check that all priors are there
    
    ## priors_names <- setdiff(names(priors), "all")
    
    ## expected_names <- grep("prior",
    ##                        names(create_config()),
    ##                        value = TRUE)
    ## expected_names <- sub("prior_", "", expected_names)
    
    ## if (!all(expected_names %in% priors_names)) {
    ##     names_missing <- expected_names[!expected_names %in% priors_names]
    ##     msg <- paste0("The following priors are missing: ",
    ##                   paste(names_missing, collapse = ", "))
    ##     stop(msg)
    ## }

    
    ## check all priors are functions

    function_or_null <- function(x) {
        is.null(x) || is.function(x)
    }
    
    is_ok <- vapply(priors, function_or_null, logical(1))
    
    if (!all(is_ok)) {
        culprits <- priors_names[!is_ok]
        msg <- paste0("The following priors are not functions: ",
                      paste(culprits, collapse = ", "))
        stop(msg)
    }

    
    ## check they all have a single parameter

    with_one_param <- function(x) {
        if(is.function(x)) {
            return (length(formalArgs(x)) == 1L)
        }
        
        return(TRUE)
    }
    
    one_arg <- vapply(priors, with_one_param, logical(1))

    if (!all(one_arg)) {
        culprits <- priors_names[!one_arg]
        msg <- paste0("The following priors dont' have a single argument: ",
                      paste(culprits, collapse = ", "))
        stop(msg)
    }
    

    ## ## Function summing all priors; we only add this one if it is missing

    ## if (is.null(priors$all)) {
    ##     priors$all <- function(param) {
    ##         sum(vapply(priors[priors_names],
    ##                    function(f) f(param), numeric(1)),
    ##             na.rm = TRUE)
    ##     }
    ## }

    priors$all <- function(param) {
        cpp_prior_all(param, config, priors$mu, priors$pi)
    }
    

    class(priors) <- c("outbreaker_priors", "list")
    return(priors)
}
