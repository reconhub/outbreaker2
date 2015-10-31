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
outbreaker.mcmc.init <- function(data, config) {
  ## CREATE EMPTY OUTPUT VECTORS ##
    size <- round(config$n.iter/config$sample.every) + 1 # +1 because initial state is always there
    post <- prior <- like <- mu <- double(size)
    ances <- as.list(integer(size))
    t.inf <- as.list(integer(size))

    ## SET CURRENT VALUES AND COUNTER ##
    current.mu <- mu[1] <- config$init.mu
    current.ances <- ances[[1]] <- config$ances
    current.t.inf <- t.inf[[1]] <- data$dates - which.max(data$f.dens) + 1
    counter <- 1L


    ## SHAPE CHAIN ##
    out <- list(size=size, post=post, like=like, prior=prior,
                ances=ances, t.inf=t.inf, mu=mu,
                current.ances=current.ances, current.t.inf=current.t.inf, current.mu=current.mu,
                counter=counter)

    ## COMPUTE INITIAL LIKE/PRIOR/POST ##
    out$like[1] <- ll.all(data=data, chain=out)
    out$prior[1] <- prior.all(chain=out)
    out$post[1] <- out$like[1] + out$prior[1]

    return(out)
} # end outbreaker.mcmc.init

