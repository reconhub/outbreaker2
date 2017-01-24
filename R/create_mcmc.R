## #' Initianilizes outputs for outbreaker
## #'
## #' This function creates initial outputs and parameter states for outbreaker.
## #'
## #' @author Thibaut Jombart (\email{t_jombart@@imperial_ac_uk})
## #'
## #' @param data a list of data items as returned by \code{outbreaker_data}
## #'
## #' @param config a list of settings as returned by \code{create_config}
## #'
## #' @export
## #'
create_mcmc <- function(data, config) {
    ## CREATE EMPTY OUTPUT VECTORS ##
    size <- round(config$n_iter/config$sample_every)
    step <- integer(size)
    post <- prior <- like <- mu <- pi <- double(size)
    alpha <- as.list(integer(size))
    t_inf <- as.list(integer(size))
    kappa <- as.list(integer(size))

    ## SET CURRENT VALUES AND COUNTER ##
    step[1] <- 1L
    current_mu <- mu[1] <- config$init_mu
    current_alpha <- alpha[[1]] <- config$init_alpha
    current_kappa <- kappa[[1]] <- config$init_kappa
    current_pi <- pi[1] <- config$init_pi
    if (is.null(config$init_t_inf)) {
        current_t_inf <- t_inf[[1]] <- data$dates - which.max(data$f_dens) + 1L
    } else {
        current_t_inf <- config$init_t_inf
    }
    counter <- 1L


    store <- list(
                size = size, step = step,
                post = post, like = like, prior = prior,
                alpha = alpha, t_inf = t_inf, mu = mu, kappa = kappa, pi = pi,
                counter = counter
                )

    current  <- list(
        alpha = current_alpha, t_inf = current_t_inf, mu = current_mu,
        kappa = current_kappa, pi = current_pi
    )
    class(current) <- c("outbreaker_param", "list")
    
    
    
    ## SHAPE CHAIN ##
    out <- list(store = store,
                current = current)
    return(out)
}

