#' outbreaker2: main function for reconstructing disease outbreaks
#'
#' The function \code{outbreaker} is the main function of the package. It runs
#' processes various inputs (data, configuration settings, custom priors,
#' likelihoods and movement functions) and explores the space of plausible
#' transmission trees of a densely sampled outbreaks.\cr
#'
#' The emphasis of 'outbreaker2' is on modularity, which enables customisation
#' of priors, likelihoods and even movements of parameters and augmented data by
#' the user. This the dedicated vignette on this topic
#' \code{vignette("outbreaker2_custom")}.\cr
#'
#' 
#'
#' @export
#'
#' @aliases outbreaker
#'
#' @rdname outbreaker
#'
#' @author Thibaut Jombart (\email{thibautjombart@@gmail.com})
#'
#' @references Jombart T, Cori A, Didelot X, Cauchemez S, Fraser C and Ferguson
#'     N (2014).  Bayesian reconstruction of disease outbreaks by combining
#'     epidemiologic and genomic data. PLoS Computational Biology.
#'
#' @seealso \code{\link{outbreaker_data}} to process input data, and
#'     \code{\link{create_config}} to process/set up parameters
#'
#' @param data a list of named items containing input data as returned by
#'     \code{\link{outbreaker_data}}
#'
#' @param config a set of settings as returned by \code{\link{create_config}}
#'
#' @param likelihoods a set of log-likelihood functions as returned by
#'     \code{\link{custom_likelihoods}}
#'
#' @param priors a set of log-prior functions as returned by
#'     \code{\link{custom_priors}}
#'
#' @param moves a set of movement functions as returned by
#'     \code{\link{custom_moves}}

#' @seealso
#'
#' \itemize{
#'
#' \item \code{\link{outbreaker_data}}: function to process input data
#'
#' \item \code{\link{create_config}}: function to create default and customise
#' configuration settings
#'
#' \item \code{\link{custom_priors}}: function to specify customised prior
#' functions
#'
#' \item \code{\link{custom_likelihoods}}: function to specify customised likelihoods
#' functions
#'
#' \item \code{\link{custom_moves}}: function to create default and customise movement
#' functions
#'
#' }
#'
#' @examples
#'
#' ## get data
#' data(fake_outbreak)
#' dat <- fake_outbreak
#'
#' \dontrun{
#' ## run outbreaker
#' out <- outbreaker(data = list(dna = dat$dna, dates = dat$onset, w_dens = dat$w),
#' config = list(n_iter = 2e4, sample_every = 200))
#' plot(out)
#' as.data.frame(out)
#'
#' ## run outbreaker, no DNA sequences
#' out2 <- outbreaker(data = list(dates = dat$onset, w_dens = w),
#' config = list(n_iter = 2e4, sample_every = 200))
#' plot(out2)
#' as.data.frame(out2)
#'
#' }
outbreaker <- function(data = outbreaker_data(),
                       config = create_config(),
                       priors = custom_priors(),
                       likelihoods = custom_likelihoods(),
                       moves = custom_moves()
                       ) {

    ## CHECKS / PROCESS DATA ##
    data <- outbreaker_data(data = data)
    
    ## CHECK / PROCESS CONFIG ##
    config <- create_config(config, data = data)

    ## ADD CONVOLUTIONS TO DATA ##
    data <- add_convolutions(data = data, config = config)

    ## PROCESS CUSTOM FUNCTIONS FOR PRIORS AND LIKELIHOOD ##
    priors <- custom_priors(priors)
    loglike <- custom_likelihoods(likelihoods)


    ## CREATE AND INITIALIZE MCMC CHAIN ##
    temp <- create_param(data = data, config = config)
    param_store <- temp$store
    param_current <- temp$current
    param_store <- outbreaker_init_mcmc(data, param_current, param_store,
                                        loglike, priors, config)


    ## here we create a list of function for moving parameters
    moves <- bind_moves(moves = moves,
                        config = config,
                        data = data,
                        likelihoods = loglike,
                        priors = priors)


    ## IMPORTS

    ## preliminary run to detect imported cases this relies on a shorter run of
    ## the MCMC, then computing the average 'global influence' (-loglike) of
    ## each data point, identifying outliers (based on fixed threshold) and
    ## marking outliers down as 'imported cases'.

    temp <- outbreaker_find_imports(moves = moves,
                                    data = data,
                                    param_current = param_current,
                                    param_store = param_store,
                                    config = config,
                                    likelihoods = loglike)
    param_current <- temp$param_current
    param_store <- temp$param_store
    config <- temp$config


    ## PERFORM MCMC

    ## procedure is the same as before, with some cases fixed as 'imported'

    param_store <- outbreaker_move(moves = moves,
                                   data = data,
                                   param_current = param_current,
                                   param_store = param_store,
                                   config = config,
                                   likelihoods = loglike,
                                   priors = priors)


    ## SHAPE RESULTS

    ## this takes the chains generated by 'outbreaker_move', stored as a list,
    ## and puts everything back together as a single data.frame; augmented data
    ## stored as vectors (e_g. 'alpha') become numbered columns of the
    ## data.frame (e_g. 'alpha_1', 'alpha_2' etc.)

    out <- outbreaker_mcmc_shape(param_store, data)

    return(out)
}

