#' Shape MCMC samples into outputs for outbreaker
#'
#' This function shapes MCMC samples for outbreaker into a data.frame with the class \code{outbreaker.chains}.
#'
#' @author Thibaut Jombart (\email{thibautjombart@@gmail.com})
#'
#' @param param a list of output items as returned by \code{outbreaker.create.mcmc}
#'
#' @param data a list of data items as returned by \code{outbreaker.data}
#'
#' @export
#'
outbreaker.mcmc.shape <- function(param, data) {
    ## unfold ancestries ##
    if (!all(vapply(param$alpha, length, integer(1))==data$N)) {
        stop("some ancestries are missing in the param")
    }
    param$alpha <- matrix(unlist(param$alpha), ncol=data$N, byrow=TRUE)
    colnames(param$alpha) <- paste("alpha", seq_len(data$N), sep=".")

    ## unfold infection dates ##
    if (!all(vapply(param$t.inf, length, integer(1))==data$N)) {
        stop("some infection dates are missing in the param")
    }
    param$t.inf <- matrix(unlist(param$t.inf), ncol=data$N, byrow=TRUE)
    colnames(param$t.inf) <- paste("t.inf", seq_len(data$N), sep=".")

    ## unfold number of generations ##
    if (!all(vapply(param$kappa, length, integer(1))==data$N)) {
        stop("some ancestries are missing in the param")
    }
    param$kappa <- matrix(unlist(param$kappa), ncol=data$N, byrow=TRUE)
    colnames(param$kappa) <- paste("kappa", seq_len(data$N), sep=".")

    ## shape data.frame and convert ##
    param <- data.frame(step=param$step,
                        post=param$post, like=param$like, prior=param$prior,
                        mu=param$mu, pi=param$pi, eps=param$eps,
                        param$alpha, param$t.inf, param$kappa)

    ## output is a data.frame containing all parameters and augmented data, with a dedicated
    ## class (for summary, plotting, etc.)
    class(param) <- c("outbreaker.chains","data.frame")
    return(param)
}
