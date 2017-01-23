

## This function creates a named list of movement functions taking a single
## argument 'param'; all the rest (e.g. likelihood, prior, posterior functions,
## config, etc) is enclosed in the functions.

bind_moves <- function(moves, config, data, likelihoods, priors) {

    out <- moves
    
   
    ## Binding:

    ## Each function needs to go through binding separately, as custom functions
    ## for likelihoods and priors often correspond to the same argument in
    ## different functions.
    
 
    ## remove move$mu if disabled
    if (!any(config$move_mu)) {
        out$mu <- NULL
    } else {
        out$mu <- bind_to_function(out$mu,
                                   data = data,
                                   config = config,
                                   custom_ll = likelihoods$genetic,
                                   custom_prior = priors$mu
                                   )
    }
    

    ## remove move$pi if disabled
    if (!any(config$move_pi)) {
        out$pi <- NULL
    } else {
        out$pi <- bind_to_function(out$pi,
                                   data = data,
                                   config = config,
                                   custom_ll = likelihoods$reporting,
                                   custom_prior = priors$pi
                                   )
    }
    

    ## remove move$alpha if no ancestry can be moved
    if (!any(config$move_alpha)) {
        out$alpha <- NULL
    } else {
        out$alpha <- bind_to_function(out$alpha,
                                      data = data,
                                      list_custom_ll = likelihoods
                                      )
    }
    

    ## remove move$t_inf if disabled
    if (!any(config$move_t_inf)) {
        out$t_inf <- NULL
    } else {
        out$t_inf <- bind_to_function(out$t_inf,
                                      data = data,
                                      list_custom_ll = likelihoods
                                      )
    }


    ## remove swap if disabled
    if (!any(config$move_swap_cases)) {
        out$swap_cases <- NULL
    } else {
        out$swap_cases <- bind_to_function(out$swap_cases,
                                      data = data,
                                      list_custom_ll = likelihoods
                                      )
    }


    ## remove move$kappa if disabled
    if (!any(config$move_kappa)) {
        out$kappa <- NULL
    } else {
        out$kappa <- bind_to_function(out$kappa,
                                      data = data,
                                      config = config,
                                      list_custom_ll = likelihoods
                                      )
    }


    ## the output is a list of movement functions with enclosed objects ##
    return(out)

}

