#' Shape MCMC samples into outputs for outbreaker
#'
#' This function shapes MCMC samples for outbreaker into a mcmc object from the package \code{coda}.
#'
#' @author Thibaut Jombart (\email{t.jombart@@imperial.ac.uk})
#'
#' @param param a list of output items as returned by \code{outbreaker.mcmc.init}
#' @param data a list of data items as returned by \code{outbreaker.data}
#'
#' @export
#'
outbreaker.mcmc.shape <- function(param, data)
{
    ## UNFOLD ANCESTRIES ##
    if(!all(sapply(param$ances, length)==data$N)) stop("some ancestries are missing in the param")
    param$ances <- matrix(unlist(param$ances), ncol=data$N, byrow=TRUE)
    colnames(param$ances) <- paste("alpha",1:data$N, sep=".")

    ## UNFOLD INFECTION DATES ##
    if(!all(sapply(param$t.inf, length)==data$N)) stop("some infection dates are missing in the param")
    param$t.inf <- matrix(unlist(param$t.inf), ncol=data$N, byrow=TRUE)
    colnames(param$t.inf) <- paste("t.inf",1:data$N, sep=".")

    ## SHAPE DATA.FRAME AND CONVERT ##
    param <- data.frame(step=param$step,
                        post=param$post, like=param$like, prior=param$prior,
                        mu=param$mu, param$ances, param$t.inf)

    ## RETURN ##
    class(param) <- c("outbreaker.chains","data.frame")
    return(param)
} # end store.mcmc
