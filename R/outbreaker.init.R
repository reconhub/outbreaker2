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
outbreaker.init <- function(data, config) {
  ## create output templates ##
    out.size <- round(config$n.iter/config$sample.every) + 1 # +1 because initial state is always there
    out.post <- out.prior <- out.like <- out.mu <- double(out.size)
    out.ances <- as.list(integer(out.size))
    out.t.inf <- as.list(integer(out.size))

    ## initialize algorithm and outputs ##
    current.mu <- out.mu[1] <- config$init.mu
    current.ances <- out.ances[[1]] <- config$ances
    current.t.inf <- out.t.inf[[1]] <- data$dates - which.max(data$f.dens) + 1
    out.like[1] <- ll.all(t.inf=out.t.inf[[1]], ances=out.ances[[1]],
                          log.w=data$log.w.dens, log.f=data$log.f.dens,
                          sampling.times=data$dates,
                          D=data$D, mu=config$init.mu, gen.length=data$L)
    out.prior[1] <- prior.all(config$init.mu)
    out.post[1] <- out.like[1] + out.prior[1]
    OUT.COUNTER <- 1L

    ## initialize pre-drawn random arrays ##
    RAND.MU <- rnorm(config$n.iter, mean=0, sd=config$sd.mu)
    RAND.ACC.MU <- log(runif(config$n.iter))
    RAND.ACC.T.INF <- log(runif(config$n.iter))
    RAND.ACC.ANCES <- log(runif(config$n.iter))
} # end outbreaker.init





