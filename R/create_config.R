#' Set and check parameter settings of outbreaker
#'
#' This function defines settings for outbreaker.  It takes a list of named
#' items as input, performs various checks, set defaults where arguments are
#' missing, and return a correct list of settings. If no input is given, it
#' returns the default settings.
#'
#' Acceptables arguments for ... are:
#'
#' \describe{
#'
#' \item{init_tree}{the tree used to initialize the MCMC. Can be either a
#' character string indicating how this tree should be computed, or a vector of
#' integers corresponding to the tree itself, where the i-th value corresponds
#' to the index of the ancestor of 'i' (i.e., \code{init.tree[i]} is the
#' ancestor of case \code{i}). Accepted character strings are "seqTrack" (uses
#' seqTrack algorithm to generate the initial tree - see function
#' \code{seqTrack} in the package \code{adegenet}), "random" (ancestor randomly
#' selected from preceding cases), and "star" (all cases coalesce to the first
#' case).  Note that for SeqTrack, all cases should have been sequenced.}
#'
#' \item{init_alpha}{a vector of integers indicating the initial values of
#' alpha, where the i-th value indicates the ancestor of case 'i'; defaults to
#' \code{NULL}, in which ancestries are defined from \code{init_tree}.}
#'
#' \item{init_kappa}{a (recycled) vector of integers indicating the initial
#' values of kappa; defaults to 1.}
#'
#' \item{init_t_inf}{a vector of integers indicating the initial values of
#' \code{t_inf}, i.e. dates of infection; defaults to \code{NULL}, in which case
#' the most likely \code{t_inf} will be determined from the delay to
#' reporting/symptoms distribution, and the dates of reporting/symptoms,
#' provided in \code{data}.}
#'
#' \item{init_mu}{initial value for the mutation rates.}
#'
#' \item{init_pi}{initial value for the reporting probability.}
#'
#' \item{init_eps}{initial value for the contact reporting coverage.}
#'
#' \item{init_tau}{initial value for the mobility paramter of unobserved cases}
#'
#' \item{init_lambda}{initial value for the non-infectious contact rate}
#'
#' \item{init_eta}{initial value for the contact sensitivity (e.g. the
#' proportion of transmission that have this type of contact)}
#'
#' \item{n_iter}{an integer indicating the number of iterations in the MCMC,
#' including the burnin period.}
#'
#' \item{move_alpha}{a vector of logicals indicating, for each case, if the
#' ancestry should be estimated ('moved' in the MCMC), or not, defaulting to
#' TRUE; the vector is recycled if needed.}
#'
#' \item{move_t_inf}{a vector of logicals indicating, for each case, if the
#' dates of infection should be estimated ('moved' in the MCMC), or not,
#' defaulting to TRUE; the vector is recycled if needed.}
#'
#' \item{move_mu}{a logical indicating whether the mutation rates
#' should be estimated ('moved' in the MCMC), or not, all defaulting to TRUE.}
#'
#' \item{move_pi}{a logical indicating whether the reporting probability
#' should be estimated ('moved' in the MCMC), or not, all defaulting to TRUE.}
#'
#' \item{move_eps}{a logical indicating whether the contact reporting coverage
#' should be estimated ('moved' in the MCMC), or not at all, defaulting to
#' TRUE.}
#'
#' \item{move_eta}{a logical indicating whether the contact sensitivity
#' should be estimated ('moved' in the MCMC), or not at all, defaulting to
#' TRUE.}
#'
#' \item{move_lambda}{a logical indicating whether the non-infectious contact
#' rate should be estimated ('moved' in the MCMC), or not at all, defaulting to
#' TRUE.}
#'
#' \item{move_kappa}{a logical indicating whether the number of generations
#' between two successive cases should be estimated ('moved' in the MCMC), or
#' not, all defaulting to TRUE.}
#'
#' \item{sample_every}{the frequency at which MCMC samples are retained for the
#' output.}
#'
#' \item{sd_mu}{the standard deviation for the Normal proposal for the mutation
#' rates.}
#'
#' \item{sd_pi}{the standard deviation for the Normal proposal for the reporting
#' probability.}
#'
#' \item{sd_eps}{the standard deviation for the Normal proposal for the
#' contact reporting coverage.}
#'
#' \item{sd_eta}{the standard deviation for the Normal proposal for the
#' contact sensitivity}
#'
#' \item{sd_lambda}{the standard deviation for the Normal proposal for the
#' non-infectious contact rate.}
#'
#' \item{prop_alpha_move}{the proportion of ancestries to move at each iteration
#' of the MCMC.}
#'
#' \item{prop_eps_move}{the proportion of iterations at which eps is moved}
#'
#' \item{prop_tau_move}{the proportion of iterations at which tau is moved}
#'
#' \item{prop_t_inf_move}{the proportion of infection dates to move at each
#' iteration of the MCMC.}
#'
#' \item{batch_size}{the size of the batch of random number pre-generated.}
#'
#' \item{paranoid}{a logical indicating if the paranoid mode should be used;
#' this mode is used for performing additional tests during outbreaker; it makes
#' computations substantially slower and is mostly used for debugging purposes.}
#'
#' \item{min_date}{earliest infection date possible, expressed as days since the
#' first sampling.}
#'
#' \item{max_kappa}{an integer indicating the largest number of generations
#' between any two linked cases; defaults to 5.}
#'
#' \item{prior_mu}{a numeric value indicating the rate of the exponential prior
#' for the mutation rate 'mu'.}
#'
#' \item{prior_pi}{a numeric vector of length 2 indicating the first and second
#' parameter of the beta prior for the reporting probability 'pi'.}
#'
#' \item{prior_eps}{a numeric vector of length 2 indicating the first and second
#' parameter of the beta prior for the contact reporting coverage 'eps'.}
#'
#' \item{prior_eta}{a numeric vector of length 2 indicating the first and second
#' parameter of the beta prior for the contact sensitivity 'eta'}
#'
#' \item{prior_lambda}{a numeric vector of length 2 indicating the first and
#' second parameter of the beta prior for the non-infectious contact rate
#' 'lambda'.}
#'
#' \item{ctd_directed}{a logical indicating if the contact tracing data is
#' directed or not. If yes, the first column represents the infector and the
#' second column the infectee. If ctd is provided as an epicontacts objects,
#' directionality will be taken from there.}
#'
#' \item{negative_si}{a logical indicating whether negative serial
#' intervals are epidemiologically possible. If not, ancestries with
#' negative serial intervals are discarded.}
#'
#' \item{genetic_model}{a character string indicating which genetic
#' model to use: one of "default", "mrca" or "phylogeny".}
#'
#' \item{pb}{a logical indicating if a progress bar should be displayed.}
#'
#' \item{init_t_onw}{initial value for the date of onward transmission;
#' defaults to \code{NULL}.}
#'
#' \item{init_place}{initial value for the place of infection;
#' defaults to \code{NULL}.}
#'
#' \item{move_swap_cases}{a logical indicating whether cases should be
#' swapped in the MCMC; defaults to TRUE.}
#'
#' \item{move_joint}{a logical indicating whether joint moves should be
#' used in the MCMC; defaults to TRUE if timed contact data are provided,
#' FALSE otherwise.}
#'
#' \item{move_model}{a logical indicating whether model moves should be
#' used in the MCMC; defaults to TRUE if timed contact data are provided,
#' FALSE otherwise.}
#'
#' \item{move_tau}{a logical indicating whether the mobility parameter
#' should be estimated ('moved' in the MCMC); defaults to TRUE if
#' timed contact data are provided, FALSE otherwise.}
#'
#' \item{swap_place}{a logical indicating whether places should be
#' swapped in the MCMC; defaults to TRUE.}
#'
#' \item{sd_tau}{the standard deviation for the Normal proposal for the
#' mobility parameter; defaults to 0.1.}
#'
#' \item{sd_t_inf}{the standard deviation for the Normal proposal for
#' the infection dates; defaults to 1.}
#'
#' \item{sd_t_onw}{the standard deviation for the Normal proposal for
#' the onward transmission dates; defaults to 5.}
#'
#' \item{prop_model_move}{the proportion of iterations at which the model
#' is moved; defaults to 0.1.}
#'
#' \item{find_import}{a logical indicating whether imported cases should
#' be identified; defaults to TRUE.}
#'
#' \item{outlier_threshold}{a numeric threshold used as a multiplier of
#' the mean influence to identify outliers (imported cases); defaults
#' to 5.}
#'
#' \item{n_iter_import}{the number of MCMC iterations used for import
#' detection; defaults to 5000.}
#'
#' \item{sample_every_import}{the sampling frequency for the import
#' detection step; defaults to 50.}
#'
#' \item{p_wrong}{the probability of error; defaults to 0.}
#'
#' \item{prior_tau}{a numeric vector of length 2 indicating the first and
#' second parameter of the beta prior for the mobility parameter 'tau'.}
#'
#' \item{return_ids}{a logical indicating whether to return case IDs;
#' defaults to FALSE.}
#'
#' }
#'
#' @param data an optional list of data items as returned by
#'     \code{outbreaker_data}; if provided, this allows for further checks of
#'     the outbreaker settings.
#'
#' @seealso \code{\link{outbreaker_data}} to check and process data for outbreaker.
#'
#' @author Thibaut Jombart (\email{thibautjombart@@gmail.com}).
#'
#' @export
#'
#' @examples
#' ## see default settings
#' create_config()
#'
#' ## change defaults
#' create_config(move_alpha = FALSE, n_iter = 2e5, sample_every = 1000)
#'
create_config <- function(..., data = NULL) {
  ## This function returns a list of configuration settings of the class
  ## 'outbreaker_config'. Arguments are passed through '...' as a list. If the
  ## list contains a single item which is already an outbreaker_config object,
  ## this one is still processed. This means all the checks are repeated, and
  ## the config is matched against data if data are provided. This should in
  ## principle allow using the same config object for several datasets. It
  ## also implicitely serves as a checking procedure for existing configs.

  config <- list(...)
  if (length(config) == 1L && is.list(config[[1]])) {
    config <- config[[1]]
  }

  have_ctd <- !is.null(data$ctd_timed)

  ## SET DEFAULTS
  defaults <- list(
    init_tree = c("seqTrack", "star", "random"),
    init_mu = 1e-4,
    init_alpha = NULL,
    init_kappa = 1,
    init_t_inf = NULL,
    init_t_onw = NULL,
    init_place = NULL,
    init_pi = 0.9,
    init_eps = 0.5,
    init_tau = 0.5,
    init_eta = 0.9,
    init_lambda = 0.05,
    move_alpha = TRUE,
    move_swap_cases = TRUE,
    move_t_inf = TRUE,
    move_joint = TRUE,
    move_model = TRUE,
    move_mu = TRUE,
    move_kappa = TRUE,
    move_pi = TRUE,
    move_tau = TRUE,
    move_eps = TRUE,
    move_eta = TRUE,
    move_lambda = TRUE,
    swap_place = TRUE,
    n_iter = 1e4,
    sample_every = 50,
    sd_mu = 0.0001,
    sd_pi = 0.1,
    sd_eps = 0.1,
    sd_eta = 0.1,
    sd_lambda = 0.05,
    sd_tau = 0.1,
    sd_t_inf = 1.0,
    sd_t_onw = 5.0,
    prop_alpha_move = 1 / 4,
    prop_eps_move = 0.05,
    prop_tau_move = 0.05,
    prop_t_inf_move = 0.2,
    prop_model_move = 0.1,
    paranoid = FALSE,
    min_date = -10,
    max_kappa = 5,
    find_import = TRUE,
    outlier_threshold = 5,
    n_iter_import = 5000,
    sample_every_import = 50,
    p_wrong = 0,
    prior_mu = c(1),
    prior_pi = c(10, 1),
    prior_tau = c(2, 2),
    prior_eps = c(1, 1),
    prior_eta = c(1, 1),
    prior_lambda = c(1, 1),
    ctd_directed = FALSE,
    negative_si = TRUE,
    genetic_model = "default",
    pb = TRUE,
    return_ids = FALSE
  )

  ## MODIFY CONFIG WITH ARGUMENTS ##
  config <- modify_defaults(defaults, config)

  ## CHECK CONFIG ##
  ## check init_tree
  if (is.character(config$init_tree)) {
    config$init_tree <- match.arg(
      config$init_tree,
      c("seqTrack", "star", "random")
    )
  }
  if (is.numeric(config$init_tree)) {
    config$init_alpha <- as.integer(config$init_tree)
  }

  ## check / process init_t_inf
  if (!is.null(config$init_t_inf)) {
    if (inherits(config$init_t_inf, "Date")) {
      config$init_t_inf <- config$init_t_inf - min(config$init_t_inf)
    }
    if (inherits(config$init_t_inf, "POSIXct")) {
      config$init_t_inf <- difftime(
        config$init_t_inf,
        min(config$init_t_inf),
        units = "days"
      )
    }
    config$init_t_inf <- as.integer(round(config$init_t_inf))
  }

  ## check init_mu
  if (!is.numeric(config$init_mu)) {
    stop("init_mu is not a numeric value")
  }
  if (config$init_mu < 0) {
    stop("init_mu is negative")
  }
  if (config$init_mu > 1) {
    stop("init_mu is greater than 1")
  }
  if (!is.finite(config$init_mu)) {
    stop("init_mu is infinite or NA")
  }

  ## check init_kappa
  if (!is.numeric(config$init_kappa)) {
    stop("init_kappa is not a numeric value")
  }
  config$init_kappa <- as.integer(round(config$init_kappa))
  if (any(config$init_kappa < 1, na.rm = TRUE)) {
    stop("init_kappa has values smaller than 1")
  }
  if (any(config$init_kappa > config$max_kappa, na.rm = TRUE)) {
    config$init_kappa[config$init_kappa > config$max_kappa] <- config$max_kappa
    warning(
      "values of init_kappa greater than max_kappa have been set to max_kappa"
    )
  }

  ## check / process init_place
  if (!is.null(config$init_place)) {
    config$init_place[] <- as.integer(config$init_place)
  }

  ## check init_pi
  if (!is.numeric(config$init_pi)) {
    stop("init_pi is not a numeric value")
  }
  if (config$init_pi < 0) {
    stop("init_pi is negative")
  }
  if (config$init_pi > 1) {
    stop("init_pi is greater than 1")
  }
  if (!is.finite(config$init_pi)) {
    stop("init_pi is infinite or NA")
  }

  ## define the total number of contact types
  n_timed <- length(data$ctd_timed_matrix)
  n_untimed <- length(data$ctd_matrix)
  n_both <- n_timed + n_untimed

  ## check init_eps
  if (any(!is.numeric(config$init_eps))) {
    stop("init_eps is not a numeric value")
  }
  if (any(config$init_eps < 0)) {
    stop("init_eps is negative")
  }
  if (any(config$init_eps > 1)) {
    stop("init_eps is greater than 1")
  }
  if (any(!is.finite(config$init_eps))) {
    stop("init_eps is infinite or NA")
  }
  if (n_both != 0 && !length(config$init_eps) %in% c(1, n_both)) {
    stop(sprintf("init_eps must be of length 1 or %d", n_both))
  }
  if (n_both != 0 && length(config$init_eps) == 1 && n_both > 1) {
    config$init_eps <- rep(config$init_eps[1], n_both)
  }

  ## check init_eta
  if (any(!is.numeric(config$init_eta))) {
    stop("init_eta is not a numeric value")
  }
  if (any(config$init_eta < 0)) {
    stop("init_eta is negative")
  }
  if (any(config$init_eta > 1)) {
    stop("init_eta is greater than 1")
  }
  if (any(!is.finite(config$init_eta))) {
    stop("init_eta is infinite or NA")
  }
  if (
    length(data$ctd_matrix) != 0 &&
      !length(config$init_eta) %in% c(1, length(data$ctd_matrix))
  ) {
    stop(sprintf("init_eta must be of length 1 or %d", length(data$ctd_matrix)))
  }
  if (
    length(data$ctd_matrix) != 0 &&
      length(config$init_eta) == 1 &&
      length(data$ctd_matrix) > 1
  ) {
    config$init_eta <- rep(config$init_eta[1], length(data$ctd_matrix))
  }

  ## check init_lambda
  if (!is.numeric(config$init_lambda)) {
    stop("init_lambda is not a numeric value")
  }

  if (any(config$init_lambda < 0)) {
    stop("init_lambda is negative")
  }
  if (any(config$init_lambda > 1)) {
    stop("init_lambda is greater than 1")
  }
  if (any(!is.finite(config$init_lambda))) {
    stop("init_lambda is infinite or NA")
  }
  if (
    length(data$ctd_matrix) != 0 &&
      !length(config$init_lambda) %in% c(1, length(data$ctd_matrix))
  ) {
    stop(sprintf(
      "init_lambda must be of length 1 or %d",
      length(data$ctd_matrix)
    ))
  }
  if (
    length(data$ctd_matrix) != 0 &&
      length(config$init_lambda) == 1 &&
      length(data$ctd_matrix) > 1
  ) {
    config$init_lambda <- rep(config$init_lambda[1], length(data$ctd_matrix))
  }

  ## check init_tau
  if (!is.numeric(config$init_tau)) {
    stop("init_tau is not a numeric value")
  }
  if (any(config$init_tau < 0)) {
    stop("init_tau is negative")
  }
  if (any(config$init_tau > 1)) {
    stop("init_tau is greater than 1")
  }
  if (any(!is.finite(config$init_tau))) {
    stop("init_tau is infinite or NA")
  }
  if (
    n_timed != 0 &&
      !length(config$init_tau) %in% c(1, n_timed)
  ) {
    stop(sprintf("init_tau must be of length 1 or %d", n_timed))
  }
  if (
    n_timed != 0 &&
      length(config$init_tau) == 1 &&
      n_timed > 1
  ) {
    config$init_tau <- rep(config$init_tau[1], n_timed)
  }

  ## check move_alpha
  if (!all(is.logical(config$move_alpha))) {
    stop("move_alpha is not a logical")
  }
  if (any(is.na(config$move_alpha))) {
    stop("move_alpha is NA")
  }

  ## check move_joint
  if (!all(is.logical(config$move_joint))) {
    stop("move_joint is not a logical")
  }
  if (any(is.na(config$move_joint))) {
    stop("move_joint is NA")
  }

  ## check move_model
  if (!all(is.logical(config$move_model))) {
    stop("move_model is not a logical")
  }
  if (any(is.na(config$move_model))) {
    stop("move_model is NA")
  }

  ## check move_swap_cases
  if (!is.logical(config$move_swap_cases)) {
    stop("move_swap_cases is not a logical")
  }
  if (is.na(config$move_swap_cases)) {
    stop("move_swap_cases is NA")
  }

  ## check move_t_inf
  if (!is.logical(config$move_t_inf)) {
    stop("move_t_inf is not a logical")
  }
  if (any(is.na(config$move_t_inf))) {
    stop("move_t_inf has NAs")
  }

  ## check move_mu
  if (!is.logical(config$move_mu)) {
    stop("move_mu is not a logical")
  }
  if (is.na(config$move_mu)) {
    stop("move_mu is NA")
  }

  ## check move_kappa
  if (!is.logical(config$move_kappa)) {
    stop("move_kappa is not a logical")
  }
  if (any(is.na(config$move_kappa))) {
    stop("move_kappa has NA")
  }

  ## check move_pi
  if (!is.logical(config$move_pi)) {
    stop("move_pi is not a logical")
  }
  if (is.na(config$move_pi)) {
    stop("move_pi is NA")
  }

  ## check move_eps
  if (!is.logical(config$move_eps)) {
    stop("move_eps is not a logical")
  }

  if (any(is.na(config$move_eps))) {
    stop("move_eps is NA")
  }
  if (n_both != 0 && !length(config$move_eps) %in% c(1, n_both)) {
    stop(sprintf("move_eps must be of length 1 or %d", n_both))
  }
  if (n_both != 0 && length(config$move_eps) == 1 && n_both > 1) {
    config$move_eps <- rep(config$move_eps[1], n_both)
  }

  ## check move_tau
  if (!is.logical(config$move_tau)) {
    stop("move_tau is not a logical")
  }
  if (any(is.na(config$move_tau))) {
    stop("move_tau is NA")
  }
  if (n_timed != 0 && !length(config$move_tau) %in% c(1, n_timed)) {
    stop(sprintf("move_tau must be of length 1 or %d", n_timed))
  }
  if (n_timed != 0 && length(config$move_tau) == 1 && n_timed > 1) {
    config$move_tau <- rep(config$move_tau[1], n_timed)
  }

  ## check move_eta
  if (!is.logical(config$move_eta)) {
    stop("move_eta is not a logical")
  }
  if (any(is.na(config$move_eta))) {
    stop("move_eta is NA")
  }
  if (n_untimed != 0 && !length(config$move_eta) %in% c(1, n_untimed)) {
    stop(sprintf("move_eta must be of length 1 or %d", n_untimed))
  }
  if (n_untimed != 0 && length(config$move_eta) == 1 && n_untimed > 1) {
    config$move_eta <- rep(config$move_eta[1], n_untimed)
  }

  ## check move_lambda
  if (!is.logical(config$move_lambda)) {
    stop("move_lambda is not a logical")
  }

  if (any(is.na(config$move_lambda))) {
    stop("move_lambda is NA")
  }
  if (n_untimed != 0 && !length(config$move_lambda) %in% c(1, n_untimed)) {
    stop(sprintf("move_lambda must be of length 1 or %d", n_untimed))
  }
  if (n_untimed != 0 && length(config$move_lambda) == 1 && n_untimed > 1) {
    config$move_lambda <- rep(config$move_lambda[1], n_untimed)
  }

  ## check n_iter
  if (!is.numeric(config$n_iter)) {
    stop("n_iter is not a numeric value")
  }
  if (config$n_iter < 2) {
    stop("n_iter is smaller than 2")
  }
  if (!is.finite(config$n_iter)) {
    stop("n_iter is infinite or NA")
  }

  ## check sample_every
  if (!is.numeric(config$sample_every)) {
    stop("sample_every is not a numeric value")
  }
  if (config$sample_every < 1) {
    stop("sample_every is smaller than 1")
  }
  if (!is.finite(config$sample_every)) {
    stop("sample_every is infinite or NA")
  }
  config$sample_every <- as.integer(config$sample_every)

  ## check sd_mu
  if (!is.numeric(config$sd_mu)) {
    stop("sd_mu is not a numeric value")
  }
  if (config$sd_mu < 1e-10) {
    stop("sd_mu is close to zero or negative")
  }
  if (!is.finite(config$sd_mu)) {
    stop("sd_mu is infinite or NA")
  }

  ## check sd_pi
  if (!is.numeric(config$sd_pi)) {
    stop("sd_pi is not a numeric value")
  }
  if (config$sd_pi < 1e-10) {
    stop("sd_pi is close to zero or negative")
  }
  if (!is.finite(config$sd_pi)) {
    stop("sd_pi is infinite or NA")
  }

  ## check sd_t_inf
  if (!is.numeric(config$sd_t_inf)) {
    stop("sd_t_inf is not a numeric value")
  }
  if (config$sd_t_inf < 1e-10) {
    stop("sd_t_inf is close to zero or negative")
  }
  if (!is.finite(config$sd_t_inf)) {
    stop("sd_t_inf is infinite or NA")
  }

  ## check sd_t_onw
  if (!is.numeric(config$sd_t_onw)) {
    stop("sd_t_onw is not a numeric value")
  }
  if (config$sd_t_onw < 1e-10) {
    stop("sd_t_onw is close to zero or negative")
  }
  if (!is.finite(config$sd_t_onw)) {
    stop("sd_t_onw is infinite or NA")
  }

  ## check sd_eps
  if (!is.numeric(config$sd_eps)) {
    stop("sd_eps is not a numeric value")
  }

  if (any(config$sd_eps < 1e-10)) {
    stop("sd_eps is close to zero or negative")
  }
  if (any(!is.finite(config$sd_eps))) {
    stop("sd_eps is infinite or NA")
  }
  if (n_both != 0 && !length(config$sd_eps) %in% c(1, n_both)) {
    stop(sprintf("sd_eps must be of length 1 or %d", n_both))
  }
  if (n_both != 0 && length(config$sd_eps) == 1 && n_both > 1) {
    config$sd_eps <- rep(config$sd_eps[1], n_both)
  }

  ## check sd_tau
  if (!is.numeric(config$sd_tau)) {
    stop("sd_tau is not a numeric value")
  }
  if (any(config$sd_tau < 1e-10)) {
    stop("sd_tau is close to zero or negative")
  }
  if (any(!is.finite(config$sd_tau))) {
    stop("sd_tau is infinite or NA")
  }
  if (n_timed != 0 && !length(config$sd_tau) %in% c(1, n_timed)) {
    stop(sprintf("sd_tau must be of length 1 or %d", n_timed))
  }
  if (n_timed != 0 && length(config$sd_tau) == 1 && n_timed > 1) {
    config$sd_tau <- rep(config$sd_tau[1], n_timed)
  }

  ## check sd_eta
  if (!is.numeric(config$sd_eta)) {
    stop("sd_eta is not a numeric value")
  }
  if (any(config$sd_eta < 1e-10)) {
    stop("sd_eta is close to zero or negative")
  }
  if (any(!is.finite(config$sd_eta))) {
    stop("sd_eta is infinite or NA")
  }
  if (n_untimed != 0 && !length(config$sd_eta) %in% c(1, n_untimed)) {
    stop(sprintf("sd_eta must be of length 1 or %d", n_untimed))
  }
  if (n_untimed != 0 && length(config$sd_eta) == 1 && n_untimed > 1) {
    config$sd_eta <- rep(config$sd_eta[1], n_untimed)
  }

  ## check sd_lambda
  if (!is.numeric(config$sd_lambda)) {
    stop("sd_lambda is not a numeric value")
  }

  if (any(config$sd_lambda < 1e-10)) {
    stop("sd_lambda is close to zero or negative")
  }
  if (any(!is.finite(config$sd_lambda))) {
    stop("sd_lambda is infinite or NA")
  }
  if (n_untimed != 0 && !length(config$sd_lambda) %in% c(1, n_untimed)) {
    stop(sprintf("sd_lambda must be of length 1 or %d", n_untimed))
  }
  if (n_untimed != 0 && length(config$sd_lambda) == 1 && n_untimed > 1) {
    config$sd_lambda <- rep(config$sd_lambda[1], n_untimed)
  }

  ## check prop_alpha_move
  if (!is.numeric(config$prop_alpha_move)) {
    stop("prop_alpha_move is not a numeric value")
  }
  if (config$prop_alpha_move < 0) {
    stop("prop_alpha_move is negative")
  }
  if (config$prop_alpha_move > 1) {
    stop("prop_alpha_move is greater than one")
  }
  if (!is.finite(config$prop_alpha_move)) {
    stop("prop_alpha_move is infinite or NA")
  }

  ## check prop_eps_move
  if (!is.numeric(config$prop_eps_move)) {
    stop("prop_eps_move is not a numeric value")
  }
  if (config$prop_eps_move < 0) {
    stop("prop_eps_move is negative")
  }
  if (config$prop_eps_move > 1) {
    stop("prop_eps_move is greater than one")
  }
  if (!is.finite(config$prop_eps_move)) {
    stop("prop_eps_move is infinite or NA")
  }

  ## check prop_tau_move
  if (!is.numeric(config$prop_tau_move)) {
    stop("prop_tau_move is not a numeric value")
  }
  if (config$prop_tau_move < 0) {
    stop("prop_tau_move is negative")
  }
  if (config$prop_tau_move > 1) {
    stop("prop_tau_move is greater than one")
  }
  if (!is.finite(config$prop_tau_move)) {
    stop("prop_tau_move is infinite or NA")
  }

  ## check prop_t_inf_move
  if (!is.numeric(config$prop_t_inf_move)) {
    stop("prop_t_inf_move is not a numeric value")
  }
  if (config$prop_t_inf_move < 0) {
    stop("prop_t_inf_move is negative")
  }
  if (config$prop_t_inf_move > 1) {
    stop("prop_t_inf_move is greater than one")
  }
  if (!is.finite(config$prop_t_inf_move)) {
    stop("prop_t_inf_move is infinite or NA")
  }

  ## check prop_model_move
  if (!is.numeric(config$prop_model_move)) {
    stop("prop_model_move is not a numeric value")
  }
  if (config$prop_model_move < 0) {
    stop("prop_model_move is negative")
  }
  if (config$prop_model_move > 1) {
    stop("prop_model_move is greater than one")
  }
  if (!is.finite(config$prop_model_move)) {
    stop("prop_model_move is infinite or NA")
  }

  ## check paranoid
  if (!is.logical(config$paranoid)) {
    stop("paranoid is not logical")
  }
  if (length(config$paranoid) != 1L) {
    stop("paranoid should be a single logical value")
  }
  if (is.na(config$paranoid)) {
    stop("paranoid is NA")
  }

  ## check min_date
  if (!is.numeric(config$min_date)) {
    stop("min_date is not numeric")
  }
  if (config$min_date >= 0) {
    stop("min_date is greater or equal to 0")
  }
  if (!is.finite(config$min_date)) {
    stop("min_date is infinite or NA")
  }

  ## check find_import
  if (!is.logical(config$find_import)) {
    stop("find_import is not logical")
  }
  if (length(config$find_import) != 1L) {
    stop("find_import should be a single logical value")
  }
  if (is.na(config$find_import)) {
    stop("find_import is NA")
  }

  ## check outlier_threshold
  if (!is.numeric(config$outlier_threshold)) {
    stop("outlier_threshold is not a numeric value")
  }
  if (any(config$outlier_threshold < 1)) {
    stop("outlier_threshold has values smaller than 1")
  }
  if (!is.finite(config$outlier_threshold)) {
    stop("outlier_threshold is infinite or NA")
  }

  ## check n_iter_import
  if (!is.numeric(config$n_iter_import)) {
    stop("n_iter_import is not a numeric value")
  }
  if (config$n_iter_import < 1000) {
    stop("n_iter is smaller than 1000")
  }
  if (!is.finite(config$n_iter_import)) {
    stop("n_iter_import is infinite or NA")
  }
  config$n_iter_import <- as.integer(config$n_iter_import)

  ## check sample_every_import
  if (!is.numeric(config$sample_every_import)) {
    stop("sample_every_import is not a numeric value")
  }
  if (config$sample_every_import < 1) {
    stop("sample_every_import is smaller than 1")
  }
  if (!is.finite(config$sample_every_import)) {
    stop("sample_every_import is infinite or NA")
  }
  config$sample_every_import <- as.integer(config$sample_every_import)

  ## check prior value for mu
  if (!is.numeric(config$prior_mu)) {
    stop("prior_mu is not a numeric value")
  }

  if (any((config$prior_mu < 0))) {
    stop("prior_mu is negative (it should be a rate)")
  }
  if (any(!is.finite(config$prior_mu))) {
    stop("prior_mu is infinite or NA")
  }

  ## check prior value for pi
  if (!all(is.numeric(config$prior_pi))) {
    stop("prior_pi has non-numeric values")
  }
  if (any(config$prior_pi < 0)) {
    stop("prior_pi has negative values")
  }
  if (length(config$prior_pi) != 2L) {
    stop("prior_pi should be a vector of length 2")
  }
  if (!all(is.finite(config$prior_pi))) {
    stop("prior_pi is has values which are infinite or NA")
  }

  ## check prior value for eps
  if (is.numeric(config$prior_eps) || is.data.frame(config$prior_eps)) {
    config$prior_eps <- matrix(config$prior_eps, ncol = 2)
  }
  if (any(!apply(config$prior_eps, 2, is.numeric))) {
    stop("prior_eps has non-numeric values")
  }
  if (any(config$prior_eps < 0)) {
    stop("prior_eps has negative values")
  }
  if (ncol(config$prior_eps) != 2L) {
    stop("prior_eps should be a matrix with two columns")
  }
  if (!all(is.finite(config$prior_eps))) {
    stop("prior_eps is has values which are infinite or NA")
  }

  if (n_both != 0 && !nrow(config$prior_eps) %in% c(1, n_both)) {
    stop(sprintf("prior_eps must be of 1 or %d rows", n_both))
  }
  if (n_both != 0 && nrow(config$prior_eps) == 1 && n_both > 1) {
    config$prior_eps <- config$prior_eps[rep(1, n_both), ]
  }

  ## check prior value for tau
  if (is.numeric(config$prior_tau) || is.data.frame(config$prior_tau)) {
    config$prior_tau <- matrix(config$prior_tau, ncol = 2)
  }
  if (any(!apply(config$prior_tau, 2, is.numeric))) {
    stop("prior_tau has non-numeric values")
  }
  if (any(config$prior_tau < 0)) {
    stop("prior_tau has negative values")
  }
  if (ncol(config$prior_tau) != 2L) {
    stop("prior_tau should be a matrix with two columns")
  }
  if (!all(is.finite(config$prior_tau))) {
    stop("prior_tau is has values which are infinite or NA")
  }
  if (n_timed != 0 && !nrow(config$prior_tau) %in% c(1, n_timed)) {
    stop(sprintf("prior_tau must be of 1 or %d rows", n_timed))
  }
  if (n_timed != 0 && nrow(config$prior_tau) == 1 && n_timed > 1) {
    config$prior_tau <- config$prior_tau[rep(1, n_timed), ]
  }

  ## check prior value for eta
  if (is.numeric(config$prior_eta) || is.data.frame(config$prior_eta)) {
    config$prior_eta <- matrix(config$prior_eta, ncol = 2)
  }
  if (any(!apply(config$prior_eta, 2, is.numeric))) {
    stop("prior_eta has non-numeric values")
  }
  if (any(config$prior_eta < 0)) {
    stop("prior_eta has negative values")
  }
  if (ncol(config$prior_eta) != 2L) {
    stop("prior_eta should be a matrix with two columns")
  }
  if (!all(is.finite(config$prior_eta))) {
    stop("prior_eta is has values which are infinite or NA")
  }
  if (n_untimed != 0 && !nrow(config$prior_eta) %in% c(1, n_untimed)) {
    stop(sprintf("prior_eta must be of 1 or %d rows", n_untimed))
  }
  if (n_untimed != 0 && nrow(config$prior_eta) == 1 && n_untimed > 1) {
    config$prior_eta <- config$prior_eta[rep(1, n_untimed), ]
  }

  ## check prior value for lambda
  if (is.numeric(config$prior_lambda) || is.data.frame(config$prior_lambda)) {
    config$prior_lambda <- matrix(config$prior_lambda, ncol = 2)
  }
  if (any(!apply(config$prior_lambda, 2, is.numeric))) {
    stop("prior_lambda has non-numeric values")
  }
  if (any(config$prior_lambda < 0)) {
    stop("prior_lambda has negative values")
  }
  if (ncol(config$prior_lambda) != 2L) {
    stop("prior_lambda should be a matrix with two columns")
  }
  if (!all(is.finite(config$prior_lambda))) {
    stop("prior_lambda is has values which are infinite or NA")
  }

  if (n_untimed != 0 && !nrow(config$prior_lambda) %in% c(1, n_untimed)) {
    stop(sprintf("prior_lambda must be of 1 or %d rows", n_untimed))
  }
  if (n_untimed != 0 && nrow(config$prior_lambda) == 1 && n_untimed > 1) {
    config$prior_lambda <- config$prior_lambda[rep(1, n_untimed), ]
  }

  if (!is.logical(config$pb)) {
    stop("pb must be a logical")
  }

  ## CHECKS POSSIBLE IF DATA IS PROVIDED ##
  if (!is.null(data)) {
    ## check initial tree
    if (is.character(config$init_tree)) {
      if (
        config$init_tree == "seqTrack" &&
          is.null(data$dna)
      ) {
        msg <- paste0(
          "Can't use seqTrack initialization with missing ",
          "DNA sequences; using a star-like tree"
        )
        message(msg)
        config$init_tree <- "star"
      }

      ## check initial tree
      if (
        config$init_tree == "seqTrack" &&
          nrow(data$dna) != data$N
      ) {
        msg <- sprintf(
          paste(
            "Can't use seqTrack initialization when",
            "numbers of sequences and cases differ",
            "(%d vs %d)"
          ),
          nrow(data$dna),
          data$N
        )
        message(msg)
        config$init_tree <- "star"
      }

      ## seqTrack init
      if (config$init_tree == "seqTrack") {
        D_temp <- data$D
        ## use strictly positive serial interval for starting tree
        can_be_ances_tmp <- outer(data$dates, data$dates, FUN = "<")
        diag(can_be_ances_tmp) <- FALSE
        D_temp[!can_be_ances_tmp] <- 1e30
        config$init_alpha <- apply(D_temp, 2, which.min)
        config$init_alpha[data$dates == min(data$dates)] <- NA
        config$init_alpha <- as.integer(config$init_alpha)
      } else if (config$init_tree == "star") {
        config$init_alpha <- rep(which.min(data$dates), length(data$dates))
        config$init_alpha[data$dates == min(data$dates)] <- NA
      } else if (config$init_tree == "random") {
        config$init_alpha <- ralpha(data$dates)
      }
    } else {
      ## if ancestries are provided
      if (length(config$init_alpha) != data$N) {
        stop("inconvenient length for init_alpha")
      }
      unknownAnces <- config$init_alpha < 1 | config$init_alpha > data$N
      if (any(stats::na.omit(unknownAnces))) {
        warning("some initial ancestries refer to unknown cases (idx<1 or >N)")
        config$init_alpha[unknownAnces] <- NA
      }
    }

    ## check initial t_inf
    if (!is.null(config$init_t_inf)) {
      if (any(config$init_t_inf >= data$dates, na.rm = TRUE)) {
        msg <- paste0(
          "Initial dates of infection come after ",
          "sampling dates / dates of onset."
        )
        stop(msg)
      }
    } else {
      ## set to most likely delay if t_inf not set
      max_like_delay <- which.max(data$f_dens)
      if (!is.finite(max_like_delay)) {
        max_like_delay <- 1L
      }
      config$init_t_inf <- as.integer(data$dates - max_like_delay)
    }

    ## recycle move_alpha only when a single value is passed; otherwise use as-is
    if (length(config$move_alpha) == 1L) {
      config$move_alpha <- rep(config$move_alpha, length.out = data$N)
    } else if (length(config$move_alpha) != data$N) {
      stop("move_alpha must be of length 1 or N (number of cases)")
    }

    ## recycle move_t_inf
    config$move_t_inf <- rep(config$move_t_inf, length.out = data$N)

    ## recycle move_kappa
    config$move_kappa <- rep(config$move_kappa, length.out = data$N)

    ## recycle init_kappa
    config$init_kappa <- rep(config$init_kappa, length.out = data$N)
    config$init_kappa[is.na(config$init_alpha)] <- NA

    ## disable moves for imported cases
    config$move_alpha[is.na(config$init_alpha)] <- FALSE
    config$move_kappa[is.na(config$init_alpha)] <- FALSE

    ## check / process init_t_onw
    if (!is.null(config$init_t_onw)) {
      if (inherits(config$init_t_onw, "Date")) {
        config$init_t_onw <- config$init_t_onw - min(config$init_t_onw)
      }
      if (inherits(config$init_t_onw, "POSIXct")) {
        config$init_t_onw <- difftime(
          config$init_t_onw,
          min(config$init_t_onw),
          units = "days"
        )
      }
      config$init_t_onw <- as.integer(round(config$init_t_onw))
    } else {
      ## set t_onw to the day before init_t_inf - this ensures that t_onw comes
      ## after the infection date of the infector
      init_t_onw <- config$init_t_inf - 1
      config$init_t_onw <- as.integer(round(init_t_onw))
    }

    if (!is.null(config$init_t_onw) && !is.null(data$ctd_timed)) {
      ## check initial t_onw
      unobs <- which(config$init_kappa > 1)
      if (any(config$init_t_onw[unobs] >= data$dates[unobs], na.rm = TRUE)) {
        msg <- paste0(
          "Initial dates of onward infection come after ",
          "sampling dates / dates of onset."
        )
        stop(msg)
      }
      if (
        any(
          config$init_t_onw[unobs] <
            config$init_t_inf[config$init_alpha[unobs]],
          na.rm = TRUE
        )
      ) {
        msg <- paste0(
          "Initial dates of onward infection come before ",
          "infection date of the infector."
        )
        stop(msg)
      }
    }

    ## disable moves for mu if no DNA sequences
    if (is.null(data$D) || nrow(data$D) < 1) {
      config$move_mu <- FALSE
    }

    ## disable moves for eps and lambda if no CTD is provided
    have_ctd <- !(is.null(data$ctd)) || !(is.null(data$ctd_timed))
    have_ctd_timed <- !(is.null(data$ctd_timed) || nrow(data$ctd_timed) < 1)

    if (!have_ctd && !have_ctd_timed) {
      config$move_eps <- config$move_eta <-
        config$move_lambda <- config$move_tau <- FALSE
    } else if (!have_ctd && have_ctd_timed) {
      config$move_lambda <- FALSE
    } else if (have_ctd && !have_ctd_timed) {
      config$move_tau <- FALSE
    }
    if (!have_ctd_timed) {
      config$move_joint <- config$move_model <- FALSE
    }
  }

  ## output is a list of checked settings with a dedicated class (for
  ## printing)

  class(config) <- c("outbreaker_config", "list")
  return(config)
}


#' @rdname create_config
#'
#' @export
#'
#' @aliases print.outbreaker_config
#'
#' @param x an \code{outbreaker_config} object as returned by \code{create_config}.
#'
#' @param ... further arguments to be passed to other methods.

print.outbreaker_config <- function(x, ...) {
  cat("\n\n ///// outbreaker settings ///\n")
  cat("\nclass:", class(x))
  cat("\nnumber of items:", length(x))

  cat("\n\n/// initialisation //\n")
  to_print <- grep("init", names(x))
  print(noquote(unlist(x[to_print])))
  x <- x[-to_print]

  cat("\n/// movements //\n")
  to_print <- unlist(sapply(c("move", "^sd"), grep, names(x)))
  print(noquote(sapply(x[to_print], as.character)))
  x <- x[-to_print]

  cat("\n/// chains //\n")
  to_print <- c("n_iter", "sample_every")
  print(noquote(sapply(x[to_print], as.character)))
  x <- x[-match(to_print, names(x))]

  cat("\n/// priors //\n")
  to_print <- grep("prior", names(x))
  print(noquote(unlist(x[to_print])))
  x <- x[-to_print]

  cat("\n/// imported cases //\n")
  to_print <- unlist(sapply(c("import", "threshold"), grep, names(x)))
  print(noquote(sapply(x[to_print], as.character)))
  x <- x[-to_print]

  if (length(x) > 0) {
    cat("\n/// other settings //\n")
    print(noquote(sapply(x, as.character)))
  }

  return(invisible(NULL))
}
