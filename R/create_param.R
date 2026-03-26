#' Initializes outputs for outbreaker
#'
#' This function creates initial outputs and parameter states for outbreaker.
#'
#' @author Thibaut Jombart (\email{thibautjombart@@gmail.com}).
#'
#' @param data A list of data items as returned by \code{outbreaker_data}, or
#' arguments passed to this function.
#'
#' @param config A list of settings as returned by \code{create_config}, or
#' arguments passed to this function.
#'
#' @export
#'
#' @aliases outbreaker_store
#'
#' @aliases outbreaker_param
#'
#' @return
#'
#' A list containing two components \code{$store} and
#' \code{$current}. \code{store} is a list with the class
#' \code{outbreaker_store}, used for storing 'saved' states of the
#' MCMC. \code{current} is a list with the class \code{outbreaker_param}, used
#' for storing 'current' states of the MCMC. \cr \cr
#'
#' \code{outbreaker_store} class content:
#' \itemize{
#'
#'  \item \code{size}: The length of the list, corresponding to the number of
#' samples saved from the MCMC.
#'
#'  \item \code{step}: A vector of integers of length \code{size}, storing the
#' steps of the MCMC corresponding to the saved samples.
#'
#'  \item \code{post}: A numeric vector of length \code{size}, storing
#' log-posterior values.
#'
#'  \item \code{like}: A numeric vector of length \code{size}, storing
#' log-likelihood values.
#'
#'  \item \code{prior}: A numeric vector of length \code{size},
#' storing log-prior values.
#'
#'  \item \code{alpha}: A list of length \code{size}. Each item of the list is
#' an integer vector of length \code{data$N}, storing indices (from 1 to N) of
#' infectors for each case.
#'
#'  \item \code{t_inf}: A list of length \code{size}. Each item of the list is
#' an integer vector of length \code{data$N}, storing dates of infections for
#' each case.
#'
#'  \item \code{mu}: A numeric vector of length \code{size}, storing values of
#' the mutation rate.
#'
#'  \item \code{kappa}: A list of length \code{size}. Each item of the list is
#' an integer vector of length \code{data$N}, storing the number of generations
#' before the last sampled ancestor for each case.
#'
#'  \item \code{pi}: A numeric vector of length \code{size}, storing values of
#' the reporting probability.
#'
#'  \item \code{eps}: A numeric vector of length \code{size}, storing values of
#' the contact reporting coverage.
#'
#'  \item \code{eta}: A numeric vector of length \code{size}, storing values of
#' the contact sensitivity.
#'
#'  \item \code{lambda}: A numeric vector of length \code{size}, storing values of
#' the non-infectious contact rate.
#'
#'  \item \code{counter}: A counter used to keep track of the current iteration
#' of the MCMC (used internally).
#'
#' }
#'
#'
#' \code{outbreaker_param} class content:
#' \itemize{
#'
#'  \item \code{alpha}: An integer vector of length \code{data$N}, storing
#' indices (from 1 to N) of infectors for each case.
#'
#'  \item \code{t_inf}: An integer vector of length \code{data$N}, storing dates
#' of infections for each case.
#'
#'  \item \code{mu}: The value of the mutation rate.
#'
#'  \item \code{kappa}: An integer vector of length \code{data$N}, storing the
#' number of generations before the last sampled ancestor for each case.
#'
#'  \item \code{pi}: The value of the reporting probability.
#'
#'  \item \code{eps}: The value of the contact reporting coverage.
#'
#'  \item \code{eta}: The value of the contact sensitivity.
#'
#'  \item \code{lambda}: The value of the non-infectious contact rate.
#'
#' }
#'
#' @examples
#'
#' ## load data
#' x <- fake_outbreak
#' data <- outbreaker_data(dates = x$sample, dna = x$dna, w_dens = x$w)
#'
#' ## modify config settings
#' config <- create_config(move_alpha = FALSE, n_iter = 2e5, sample_every = 1000)
#'
#' ## create param object
#' param <- create_param(data = data, config = config)
#'
create_param <- function(data = outbreaker_data(),
                         config = create_config()) {
  ## CREATE EMPTY OUTPUT VECTORS ##
  size <- round(config$n_iter/config$sample_every)
  step <- integer(size)
  post <- prior <- like <- mu <- pi <- tau <- double(size)
  eps <- as.list(integer(size))
  tau <- as.list(integer(size))
  eta <- as.list(integer(size))
  lambda <- as.list(integer(size))
  alpha <- as.list(integer(size))
  t_inf <- as.list(integer(size))
  t_onw <- as.list(integer(size))
  kappa <- as.list(integer(size))

  ## SET CURRENT VALUES AND COUNTER ##
  step[1] <- 1L
  current_mu <- mu[1] <- config$init_mu
  current_alpha <- alpha[[1]] <- config$init_alpha
  current_kappa <- kappa[[1]] <- config$init_kappa
  current_pi <- pi[1] <- config$init_pi
  current_tau <- tau[[1]] <- config$init_tau
  current_eps <- eps[[1]] <- config$init_eps
  current_eta <- eta[[1]] <- config$init_eta
  current_lambda <- lambda[[1]] <- config$init_lambda

  # Set values for p_wrong, eps and lambda away from 0 and 1 to
  # prevent -Inf likelihood in initialisation phase
  is_equal <- function(x, val)
    isTRUE(all.equal(x, val, tolerance = 1e-16))
  eps_one <- vapply(current_eps, is_equal, TRUE, 1)
  if (any(eps_one)) current_eps[eps_one] <- 1 - 1e-15
  tau_one <- vapply(current_tau, is_equal, TRUE, 1)
  if (any(tau_one)) current_tau[tau_one] <- 1 - 1e-15

  if (is.null(config$init_t_inf)) {
    current_t_inf <- t_inf[[1]] <- data$dates - which.max(data$f_dens) + 1L
  } else {
    current_t_inf <- t_inf[[1]] <- config$init_t_inf
  }
  if (is.null(config$init_t_onw)) {
    tmp_t_onw <- round(config$init_t_inf - sum(data$w_dens*seq_along(data$w_dens))/2)
    tmp_t_onw[config$init_kappa == 1] <- -1000
    current_t_onw <- t_onw[[1]] <- as.integer(tmp_t_onw)
  } else {
    tmp_t_onw <- config$init_t_onw
    tmp_t_onw[config$init_kappa == 1] <- -1000
    current_t_onw <- t_onw[[1]] <- as.integer(tmp_t_onw)
  }

  ## Calculate initial place transition probabilities
  if (!is.null(data$ctd_timed)) {
    trans_mat <- list()
    ## number of untimed contact types for indexing
    n_u <- length(data$ctd_matrix)
    for(i in seq_along(data$p_trans)) {
      trans_mat[[i]] <- get_transition_mat(data$pp_trans_adj[[i]],
                                           data$pp_place[[i]],
                                           data$pp_place_adj[[i]],
                                           current_eps[i+n_u],
                                           current_tau[i],
                                           data$prop_place_observed[i],
                                           config$max_kappa)
    }
  } else {
    trans_mat <- list()
  }

  ## ## Calculate initial ancestor tree if genetic sequences provided
  if(!is.null(data$dna) & config$genetic_model == "mrca") {
    ancestors <- matrix(0, data$N, data$N)
    ancestors <- cpp_find_ancestors(current_alpha, ancestors, NULL)
    mrca <- t(apply(
      data$dna_combn, 1,
      function(x) cpp_find_mrca(x[1], x[2], ancestors)
    ))
  } else {
    ancestors <- mrca <- matrix(0, 0, 0)
  }

  counter <- 1L

  store <- list(
    size = size, step = step,
    post = post, like = like, prior = prior,
    alpha = alpha, t_inf = t_inf, t_onw = t_onw,
    mu = mu, kappa = kappa, pi = pi, tau = tau,
    eps = eps, eta = eta, lambda = lambda,
    counter = counter
  )

  class(store) <- c("outbreaker_store", "list")

  current  <- list(alpha = current_alpha, t_inf = current_t_inf,
                   t_onw = current_t_onw, mu = current_mu,
                   kappa = current_kappa, pi = current_pi, tau = current_tau,
                   eps = current_eps, eta = current_eta,
                   lambda = current_lambda,
                   trans_mat = trans_mat,
                   ancestors = ancestors, mrca = mrca)
  class(current) <- c("outbreaker_param", "list")


  ## SHAPE CHAIN ##
  out <- list(store = store, current = current)
  return(out)

}
