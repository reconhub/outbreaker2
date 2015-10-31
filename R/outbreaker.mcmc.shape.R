#' Shape MCMC samples into outputs for outbreaker
#'
#' This function shapes MCMC samples for outbreaker into a mcmc object from the package \code{coda}.
#'
#' @author Thibaut Jombart (\email{t.jombart@@imperial.ac.uk})
#'
#' @param chain a list of output items as returned by \code{outbreaker.mcmc.init}
#'
#' @export
#'
outbreaker.mcmc.shape <- function(chain, data) {
    ## UNFOLD ANCESTRIES ##
    chain$ances <- matrix(unlist(chain$ances), ncol=data$N, byrow=TRUE)
    colnames(chain$ances) <- paste("alpha",1:data$N, sep=".")

    ## UNFOLD INFECTION DATES ##
    chain$t.inf <- matrix(unlist(chain$t.inf), ncol=data$N, byrow=TRUE)
    colnames(chain$inf) <- paste("t.inf",1:data$N, sep=".")

    ## SHAPE DATA.FRAME AND CONVERT ##
    chain <- data.frame(post=chain$post, like=chain$like, prior=chain$prior,
                      mu=chain$mu, chain$ances, chain$t.inf)
    chain <- coda::mcmc(chain)

    ## RETURN ##
    return(chain)
} # end store.mcmc

