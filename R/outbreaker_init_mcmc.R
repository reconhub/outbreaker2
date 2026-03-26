#' Initialise MCMC storage for outbreaker
#'
#' Computes the initial log-likelihood, log-prior, and log-posterior for the
#' starting parameter state and writes them into \code{param_store}.
#'
#' @param data A list of data items as returned by \code{\link{outbreaker_data}}.
#'
#' @param param_current A list of current parameter states as returned by
#'   \code{\link{create_param}}.
#'
#' @param param_store A list of stored parameter states for MCMC output.
#'
#' @param loglike A list of log-likelihood components as returned by
#'   \code{\link{custom_likelihoods}}.
#'
#' @param priors A list of prior functions as returned by
#'   \code{\link{custom_priors}}.
#'
#' @param config A list of configuration settings as returned by
#'   \code{\link{create_config}}.
#'
#' @return The \code{param_store} list, with \code{like[1]}, \code{prior[1]},
#'   and \code{post[1]} filled.
#'
#' @author Thibaut Jombart (\email{thibautjombart@@gmail.com}).
#'
#' @seealso \code{\link{outbreaker}}, \code{\link{create_param}}.
#'
#' @export
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
