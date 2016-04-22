#' Shape MCMC samples into outputs for outbreaker
#'
#' This function shapes MCMC samples for outbreaker into a data.frame with the class \code{outbreaker.chains}.
#'
#' @author Thibaut Jombart (\email{thibautjombart@@gmail.com})
#'
#' @param param a list of output items as returned by \code{outbreaker.create.mcmc}
#' @param data a list of data items as returned by \code{outbreaker.data}
#'
#' @export
#'
outbreaker.mcmc.shape <- function(param, data){
    ## UNFOLD ANCESTRIES ##
    if (!all(sapply(param$alpha, length)==data$N)) {
        stop("some ancestries are missing in the param")
    }
    param$alpha <- matrix(unlist(param$alpha), ncol=data$N, byrow=TRUE)
    colnames(param$alpha) <- paste("alpha", seq_len(data$N), sep=".")

    ## UNFOLD INFECTION DATES ##
    if (!all(sapply(param$t.inf, length)==data$N)) {
        stop("some infection dates are missing in the param")
    }
    param$t.inf <- matrix(unlist(param$t.inf), ncol=data$N, byrow=TRUE)
    colnames(param$t.inf) <- paste("t.inf", seq_len(data$N), sep=".")

    ## UNFOLD NUMBER OF GENERATIONS ##
    if (!all(sapply(param$kappa, length)==data$N)) {
        stop("some ancestries are missing in the param")
    }
    param$kappa <- matrix(unlist(param$kappa), ncol=data$N, byrow=TRUE)
    colnames(param$kappa) <- paste("kappa", seq_len(data$N), sep=".")

    ## SHAPE DATA.FRAME AND CONVERT ##
    param <- data.frame(step=param$step,
                        post=param$post, like=param$like, prior=param$prior,
                        mu=param$mu, pi=param$pi,
                        param$alpha, param$t.inf, param$kappa)

    ## RETURN ##
    class(param) <- c("outbreaker.chains","data.frame")
    return(param)
} # end store.mcmc
