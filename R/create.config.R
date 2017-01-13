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
#' \item{init.tree}{the tree used to initialize the MCMC. Can be either a
#' character string indicating how this tree should be computed, or a vector of
#' integers corresponding to the tree itself, where the i-th value corresponds
#' to the index of the ancestor of 'i' (i.e., \code{init.tree[i]} is the
#' ancestor of case \code{i}). Accepted character strings are "seqTrack" (uses
#' seqTrack algorithm to generate the initial tree - see function
#' \code{seqTrack} in the package \code{adegenet}), "random" (ancestor randomly
#' selected from preceding cases), and "star" (all cases coalesce to the first
#' case).  Note that for SeqTrack, all cases should have been sequenced.}
#'
#' \item{n.iter}{an integer indicating the number of iterations in the MCMC,
#' including the burnin period}
#'
#' \item{init.mu}{initial value for the mutation rates}
#'
#' \item{init.kappa}{a (recycled) vector of integers indicating the initial
#' values of kappa; defaults to 1.}
#'
#' \item{init.pi}{initial value for the reporting probability}
#'
#' \item{move.alpha}{a vector of logicals indicating, for each case, if the
#' ancestry should be estimated ('moved' in the MCMC), or not, defaulting to
#' TRUE; the vector is recycled if needed.}
#'
#' \item{move.t_inf}{a vector of logicals indicating, for each case, if the
#' dates of infection should be estimated ('moved' in the MCMC), or not,
#' defaulting to TRUE; the vector is recycled if needed.}
#'
#' \item{move.mu}{a logical indicating whether the mutation rates
#' should be estimated ('moved' in the MCMC), or not, all defaulting to TRUE.}
#'
#' \item{move.pi}{a logical indicating whether the reporting probability
#' should be estimated ('moved' in the MCMC), or not, all defaulting to TRUE.}
#'
#' \item{move.kappa}{a logical indicating whether the number of generations
#' between two successive cases should be estimated ('moved' in the MCMC), or
#' not, all defaulting to TRUE.}
#'
#' \item{move_pi}{a logical indicating whether the reporting probability
#' should be estimated ('moved' in the MCMC), or not, all defaulting to TRUE.}
#'
#' \item{n_iter}{the number of iterations of the MCMC}
#'
#' \item{sample_every}{the frequency at which MCMC samples are retained for the
#' output}
#'
#' \item{sd_mu}{the standard deviation for the Normal proposal for the mutation
#' rates}
#'
#' \item{sd_pi}{the standard deviation for the Normal proposal for the reporting
#' probability}
#'
#' \item{prop_alpha_move}{the proportion of ancestries to move at each iteration
#' of the MCMC}
#'
#' \item{prop_t_inf_move}{the proportion of infection dates to move at each
#' iteration of the MCMC}

#' \item{batch_size}{the size of the batch of random number pre-generated}
#'
#' \item{paranoid}{a logical indicating if the paranoid mode should be used;
#' this mode is used for performing additional tests during outbreaker; it makes
#' computations substantially slower and is mostly used for debugging purposes.}
#'
#' \item{min_date}{earliest infection date possible, expressed as days since the
#' first sampling}
#'
#' \item{max_kappa}{an integer indicating the largest number of generations
#' between any two linked cases; defaults to 5}
#'
#' \item{prior_mu}{a numeric value indicating the rate of the exponential prior
#' for the mutation rate 'mu'}
#'
#' \item{prior_pi}{a numeric vector of length 2 indicating the first and second
#' parameter of the beta prior for the reporting probability 'pi'}
#'
#' }
#'
#' @param ... settings to be passed to outbreaker
#'
#' @param data an optional list of data items as returned by \code{outbreaker_data}; if provided,
#' this allows for further checks of the outbreaker setings.
#'
#' @param config a previous set of settings, as returned by \code{create_config}
#'
#' @seealso outbreaker_data to check and process data for outbreaker
#'
#' @author Thibaut Jombart (\email{t.jombart@@imperial.ac.uk})
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
#'
#' @importFrom utils modifyList
#'
create_config <- function(..., data = NULL) {
    
    ## Arguments are passed through '...' as a list. If the list contains a
    ## single item which is already an outbreaker_config object, this one is
    ## still processed. This means all the checks are repeated, and the config
    ## is matched against data if data are provided. This should in principle
    ## allow using the same config object for several datasets. It also
    ## implicitely serves as a checking procedure for existing configs.
    
    config <- list(...)
    if (length(config) == 1L && is.list(config[[1]])) {
        config <- config[[1]]
        ## return(config[[1]])
    }

    ## SET DEFAULTS
    defaults <- list(init_tree = c("seqTrack","star","random"),
                     init_mu = 1e-4,
                     init_t_inf = NULL,
                     init_alpha = NULL,
                     init_kappa = 1,
                     init_pi = 0.9,
                     move_alpha = TRUE, move_swap_cases = TRUE, move_t_inf = TRUE,
                     move_mu = TRUE, move_kappa = TRUE, move_pi = TRUE,
                     n_iter = 1e4, sample_every = 200, sd_mu = 0.0001, sd_pi = 0.1,
                     prop_alpha_move = 1/4,
                     prop_t_inf_move = 0.2,
                     batch_size = 1e6,
                     paranoid = FALSE,
                     min_date = -10,
                     max_kappa = 5,
                     find_import = TRUE,
                     outlier_threshold = 5,
                     n_iter_import = 5000,
                     sample_every_import = 100,
                     prior_mu = 1000,
                     prior_pi = c(10,1))

    ## MODIFY CONFIG WITH ARGUMENTS ##
    config <- modify_defaults(defaults, config)

    ## CHECK CONFIG ##
    ## check init_tree
    if (is.character(config$init_tree)) {
        config$init_tree <- match.arg(config$init_tree, c("seqTrack","star","random"))
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
            config$init_t_inf <- difftime(config$init_t_inf,
                                          min(config$init_t_inf), units="days")
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
    if (any(config$init_kappa < 1)) {
        stop("init_kappa has values smaller than 1")
    }
    if (any(config$init_kappa > config$max_kappa)) {
        config$init_kappa[config$init_kappa > config$max_kappa] <- config$max_kappa
        warning("values of init_kappa greater than max_kappa have been set to max_kappa")
    }
    if (!all(is.finite(config$init_kappa))) {
        stop("init_kappa has values which are infinite or NA")
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

    ## check move_alpha
    if (!is.logical(config$move_alpha)) {
        stop("move_alpha is not a logical")
    }
    if (is.na(config$move_alpha)) {
        stop("move_alpha is NA")
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
    if (is.na(config$move_t_inf)) {
        stop("move_t_inf is NA")
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
    if (is.na(config$move_kappa)) {
        stop("move_kappa is NA")
    }

    ## check move_pi
    if (!is.logical(config$move_pi)) {
        stop("move_pi is not a logical")
    }
    if (is.na(config$move_pi)) {
        stop("move_pi is NA")
    }

    ## check n_iter
    if (!is.numeric(config$n_iter)) {
        stop("n_iter is not a numeric value")
    }
    if (config$n_iter < 2 ) {
        stop("n_iter is smaller than 2")
    }
    if (!is.finite(config$n_iter)) {
        stop("n_iter is infinite or NA")
    }


    ## check sample_every
    if (!is.numeric(config$sample_every)) {
        stop("sample_every is not a numeric value")
    }
    if (config$sample_every < 1 ) {
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

    ## check prop_alpha_move
    if (!is.numeric(config$prop_alpha_move)) {
        stop("prop_alpha_move is not a numeric value")
    }
    if (config$prop_alpha_move < 0 ) {
        stop("prop_alpha_move is negative")
    }
    if (config$prop_alpha_move > 1 ) {
        stop("prop_alpha_move is greater than one")
    }
    if (!is.finite(config$prop_alpha_move)) {
        stop("prop_alpha_move is infinite or NA")
    }

    ## check prop_t_inf_move
    if (!is.numeric(config$prop_t_inf_move)) {
        stop("prop_t_inf_move is not a numeric value")
    }
    if (config$prop_t_inf_move < 0 ) {
        stop("prop_t_inf_move is negative")
    }
    if (config$prop_t_inf_move > 1 ) {
        stop("prop_t_inf_move is greater than one")
    }
    if (!is.finite(config$prop_t_inf_move)) {
        stop("prop_t_inf_move is infinite or NA")
    }

    ## check batch_size
    if (!is.numeric(config$batch_size)) {
        stop("batch_size is not numeric")
    }
    if (config$batch_size < 1) {
        stop("batch_size is less than 1")
    }
    if (!is.finite(config$batch_size)) {
        stop("batch_size is infinite or NA")
    }
    config$batch_size <- as.integer(config$batch_size)

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
    if (config$prior_mu < 0) {
        stop("prior_mu is negative (it should be a rate)")
    }
    if (!is.finite(config$prior_mu)) {
        stop("prior_mu is infinite or NA")
    }

    ## check prior value for pi
    if (!all(is.numeric(config$prior_pi))) {
        stop("prior_pi has non-numeric values")
    }
    if (any(config$prior_pi < 0)) {
        stop("prior_pi has negative values")
    }
    if (length(config$prior_pi)!=2L) {
        stop("prior_pi should be a vector of length 2")
    }
   if (!all(is.finite(config$prior_pi))) {
        stop("prior_pi is has values which are infinite or NA")
    }


    ## CHECKS POSSIBLE IF DATA IS PROVIDED ##
    if (!is.null(data)) {
        ## check initial tree
        if (is.character(config$init_tree)) {
            if (config$init_tree=="seqTrack" && is.null(data$dna)) {
                message("Can't use seqTrack initialization with missing DNA sequences; using a star-like tree")
                config$init_tree <- "star"
            }

            ## seqTrack init
            if (config$init_tree=="seqTrack") {
                D_temp <- data$D
                D_temp[!data$can_be_ances] <- 1e30
                config$init_alpha <- apply(D_temp,2,which.min)
                config$init_alpha[data$dates==min(data$dates)] <- NA
                config$init_alpha <- as.integer(config$init_alpha)
            } else if (config$init_tree=="star") {
                config$init_alpha <- rep(which.min(data$dates), length(data$dates))
                config$init_alpha[data$dates==min(data$dates)] <- NA
            } else if (config$init_tree=="random") {
                config$init_alpha <- ralpha(data$dates)
            }
        } else { ## if ancestries are provided
            if (length(config$init_alpha) != data$N) {
                stop("inconvenient length for init_alpha")
            }
            unknownAnces <- config$init_alpha<1 | config$init_alpha>data$N
            if (any(stats::na.omit(unknownAnces))) {
                warning("some initial ancestries refer to unknown cases (idx<1 or >N)")
                config$init_alpha[unknownAnces] <- NA
            }
        }

        ## check initial t_inf
        if (!is.null(config$init_t_inf)) {
            if (any(config$init_t_inf >= data$dates)) {
                stop("Initial dates of infection come after sampling dates / dates of onset.")
            }
        }

        ## recycle move_alpha
        config$move_alpha <- rep(config$move_alpha, length.out = data$N)

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

        ## disable moves for mu if no DNA sequences
        if(is.null(data$D) || nrow(data$D)<1) {
            config$move_mu <- FALSE
        }
    }

    ## output is a list of checked settings with a dedicated class (for printing)
    class(config) <- c("outbreaker_config", "list")
    return(config)
}

