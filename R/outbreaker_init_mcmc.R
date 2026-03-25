#' Initialises outputs for outbreaker
#'
#' This function creates initial outputs and parameter states for outbreaker.
#'
#' @author Thibaut Jombart (\email{thibautjombart@@gmail.com})
#'
#' @param data a list of data items as returned by \code{outbreaker_data}.
#'
#' @param param_current a list of current parameter states as returned by
#'   \code{create_param}.
#'
#' @param param_store a list of stored parameter states for MCMC output.
#'
#' @param loglike a list of loglikelihood functions with enclosed data
#'   as returned by \code{custom_likelihoods}.
#'
#' @param priors a list of prior functions with enclosed parameters as
#'   returned by \code{custom_priors}.
#'
#' @param config a list of configuration settings as returned by
#'   \code{create_config}.
#'
#' @export
#'
outbreaker_init_mcmc <- function(data, param_current, param_store,
                                 loglike, priors, config) {

  ## COMPUTE INITIAL LIKE/PRIOR/POST ##
  param_store$like[1] <- cpp_ll_all(data, param_current, NULL, loglike)
  param_store$prior[1] <- cpp_prior_all(param_current, config, priors)
  param_store$post[1] <- param_store$like[1] + param_store$prior[1]

  if (is.infinite(param_store$post[1])) {
    stop("Likelihood of initial parameter state is -Inf")
  }

  return(param_store)

}
