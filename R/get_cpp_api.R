
#' Access internal C++ rountines used in outbreaker2
#'
#' This function returns an environment containing all C++ functions (bound to R
#' using Rcpp) used for priors, likelihoods, and movements of parameters in
#' outbreaker2.
#'
#' @export
#'
#' @author Thibaut Jombart (\email{thibautjombart@@gmail.com})
#'
#' @return
#'
#' An environment containing Rcpp bindings to C++ functions used internally in
#' outbreaker2. All functions are named as cpp_[type]_[component], where 'type'
#' can be:
#'
#' \itemize{
#'
#' \item 'prior': prior computation
#'
#' \item 'll':  likelihood computation
#'
#' \item 'move': movement of parameters or augmented data
#'
#' }
#'
#' And where 'component' can be:
#'
#' \itemize{
#'
#' \item 'mu': (parameter) mutation rate
#'
#' \item 'pi': (parameter) reporting probability
#'
#' \item 'lambda': (parameter) non-infectious contact rate
#'
#' \item 'eps': (parameter) contact reporting coverage
#' 
#' \item 'alpha': (augmented data) ancestries of the cases
#'
#' \item 'kappa': (augmented data) generation between cases on transmission
#' chains
#'
#' \item 't_inf': (augmented data) dates of infections
#'
#' \item 'timing_infections': (likelihood component) timing between infectors and
#' infectees
#'
#' \item 'timing_sampling': (likelihood component) timing between infection and
#' reporting / isolation
#'
#' \item 'timing': (likelihood component) sum of the two timing components
#' 
#' \item 'genetic': (likelihood component) genetic diversity accumulated along
#' transmission chains
#' 
#' \item 'reporting': (likelihood component) reporting of cases
#' 
#' \item 'all': (likelihood component) sum of all likelihood components
#'
#' \item 'swap_cases': (type of movement) swap infectors and infectees on a
#' transmission tree 
#'
#' }
#'
#' For a description of the arguments of these functions, see the Rcpp_API
#' vignette (\code{vignette("Rcpp_API", package = "outbreaker2")}).
#'
#' 
#' @examples
#'
#' ## get functions in an environment
#' api <- get_cpp_api()
#' api
#'
#' ## check content of the environment
#' ls(api)
#'
#' ## test the prior for mu
#' args(api$cpp_prior_mu)
#'
#' config <- create_config()
#'
#' api$cpp_prior_mu(list(mu = 0.00123), config)
#' 
#' dexp(0.00123, rate = config$prior_mu, log = TRUE)
#' 
get_cpp_api <- function() {
    pkg_env <- asNamespace("outbreaker2")
    regxp <- "^cpp_(ll|prior|move)"
    names_cpp_functions <- sort(ls(envir = pkg_env, pattern = regxp))

    out_env <- new.env()
    
    for (e in names_cpp_functions) {
        f <- get(e, envir = pkg_env)
        assign(e, f, envir = out_env)
    }

    return(out_env)
}
