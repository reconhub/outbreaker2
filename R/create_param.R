#' Initianilizes outputs for outbreaker
#'
#' This function creates initial outputs and parameter states for outbreaker.
#'
#' @author Thibaut Jombart (\email{thibautjombart@@gmail.com})
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
#'  \item \code{size}: The length of the list, corresponding to the number of samples saved from the MCMC.
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
#'  \item \code{eps}: Finlay - complete
#'
#'  \item \code{lambda}: Finlay - complete
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
#'  \item \code{eps}: Finlay - complete
#'
#'  \item \code{lambda}: Finlay - complete
#'
#' }
#'
#' 
create_param <- function(data = outbreaker_data(),
                         config = create_config()) {
    ## CREATE EMPTY OUTPUT VECTORS ##
    size <- round(config$n_iter/config$sample_every)
    step <- integer(size)
    post <- prior <- like <- mu <- pi <- eps <- lambda <- double(size)
    alpha <- as.list(integer(size))
    t_inf <- as.list(integer(size))
    kappa <- as.list(integer(size))

    ## SET CURRENT VALUES AND COUNTER ##
    step[1] <- 1L
    current_mu <- mu[1] <- config$init_mu
    current_alpha <- alpha[[1]] <- config$init_alpha
    current_kappa <- kappa[[1]] <- config$init_kappa
    current_pi <- pi[1] <- config$init_pi
    current_eps <- eps[1] <- config$init_eps
    current_lambda <- lambda[1] <- config$init_lambda
    if (is.null(config$init_t_inf)) {
        current_t_inf <- t_inf[[1]] <- data$dates - which.max(data$f_dens) + 1L
    } else {
        current_t_inf <- config$init_t_inf
    }
    counter <- 1L


    store <- list(
                size = size, step = step,
                post = post, like = like, prior = prior,
                alpha = alpha, t_inf = t_inf, mu = mu, kappa = kappa, pi = pi,
                eps = eps, lambda = lambda,
                counter = counter
                )
    class(store) <- c("outbreaker_store", "list")

    
    current  <- list(
        alpha = current_alpha, t_inf = current_t_inf, mu = current_mu,
        kappa = current_kappa, pi = current_pi, eps = current_eps,
        lambda = current_lambda
    )
    class(current) <- c("outbreaker_param", "list")



    ## SHAPE CHAIN ##
    out <- list(store = store,
                current = current)
    return(out)
}

