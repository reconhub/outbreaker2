#' Initianilizes outputs for outbreaker
#'
#' This function creates initial outputs and parameter states for outbreaker.
#'
#' @author Thibaut Jombart (\email{t.jombart@@imperial.ac.uk})
#'
#' @param data a list of data items as returned by \code{outbreaker.data}
#'
#' @param config a list of settings as returned by \code{outbreaker.config}
#'
#' @export
#'
outbreaker.create.mcmc <- function(data, config) {
    ## CREATE EMPTY OUTPUT VECTORS ##
    size <- round(config$n.iter/config$sample.every)
    step <- integer(size)
    post <- prior <- like <- mu <- pi <- double(size)
    alpha <- as.list(integer(size))
    t.inf <- as.list(integer(size))
    kappa <- as.list(integer(size))

    ## SET CURRENT VALUES AND COUNTER ##
    step[1] <- 1L
    current.mu <- mu[1] <- config$init.mu
    current.alpha <- alpha[[1]] <- config$init.alpha
    current.kappa <- kappa[[1]] <- config$init.kappa
    current.pi <- pi[1] <- config$init.pi
    if (is.null(config$init.t.inf)) {
        current.t.inf <- t.inf[[1]] <- data$dates - which.max(data$f.dens) + 1L
    } else {
        current.t.inf <- config$init.t.inf
    }
    counter <- 1L


    ## SHAPE CHAIN ##
    out <- list(store = list(
                size = size, step = step,
                post = post, like = like, prior = prior,
                alpha = alpha, t.inf = t.inf, mu = mu, kappa = kappa, pi = pi
                ),
                current = list(
                alpha = current.alpha, t.inf = current.t.inf, mu = current.mu,
                kappa = current.kappa, pi = current.pi,
                counter = counter
                )
                )
    return(out)
}

