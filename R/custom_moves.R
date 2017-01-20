

## This function creates a named list of movement functions taking a single
## argument 'param'; all the rest (e.g. likelihood, prior, posterior functions,
## config, etc) is enclosed in the functions.

create_moves <- function(config, data, densities) {

    ## These are all the functions generating various movement functions; we
    ## list them by alphabetic order.

    ## SET DEFAULTS ##
    default_functions <- list(mu = make_move_mu,
                     t_inf = make_move_t_inf,
                     alpha = make_move_alpha,
                     swap_cases = make_move_swap_cases,
                     pi = make_move_pi,
                     kappa = make_move_kappa
                     )

    out <- lapply(default_functions, function(f) f(config, data, densities))


    ## REMOVE FUNCTIONS IF MOVEMENTS DISABLED ##
    ## remove move$alpha if no ancestry can be moved
    if (!any(config$move_alpha)) {
        out$alpha <-  out$swap_cases <- NULL
    }

    ## remove move$t_inf if disabled
    if (!any(config$move_t_inf)) {
        out$t_inf <-  out$swap_cases <- NULL
    }

    ## remove move$mu if disabled
    if (!any(config$move_mu)) {
        out$mu <- NULL
    }

    ## remove swap if disabled, or if some alpha/t_inf cannot be moved
    if (!any(config$move_swap_cases) ||
       !any(config$move_alpha) ||
       !any(config$move_t_inf)) {
        out$swap_cases <- NULL
    }

    ## remove move$pi if disabled
    if (!any(config$move_pi)) {
        out$pi <- NULL
    }

    ## remove move$kappa if disabled
    if (!any(config$move_kappa)) {
        out$kappa <- NULL
    }


    ## the output is a list of movement functions with enclosed stuff ##
    return(out)

}

