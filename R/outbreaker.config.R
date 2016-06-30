#' Set and check parameter settings of outbreaker
#'
#' This function defines settings for outbreaker.
#' It takes a list of named items as input, performs various checks, set defaults where arguments are missing, and return a correct list of settings. If no input is given, it returns the default settings.
#'
#' Acceptables arguments for ... are:
#'
#' \describe{
#'
#' \item{init.tree}{the tree used to initialize the MCMC. Can be either a character string
#' indicating how this tree should be computed, or a vector of integers corresponding to the tree
#' itself, where the i-th value corresponds to the index of the ancestor of 'i' (i.e.,
#' \code{init.tree[i]} is the ancestor of case \code{i}). Accepted character strings are "seqTrack"
#' (uses seqTrack algorithm to generate the initial tree - see function \code{seqTrack} in the
#' package \code{adegenet}), "random" (ancestor randomly selected from preceding cases), and "star"
#' (all cases coalesce to the first case).  Note that for SeqTrack, all cases should have been
#' sequenced.}
#'
#' \item{n.iter}{an integer indicating the number of iterations in the MCMC,
#' including the burnin period}
#'
#' \item{init.mu}{initial value for the mutation rates}
#'
#' \item{init.kappa}{a (recycled) vector of integers indicating the initial values of kappa; defaults to 1.}
#'
#' \item{init.pi}{initial value for the reporting probability}
#'
#' \item{init.eps}{initial value for the contact reporting probability}
#'
#' \item{move.alpha}{a vector of logicals indicating, for each case, if the ancestry should be
#' estimated ('moved' in the MCMC), or not, defaulting to TRUE; the vector is recycled if needed.}
#'
#' \item{move.t.inf}{a vector of logicals indicating, for each case, if the dates of infection
#' should be estimated ('moved' in the MCMC), or not, defaulting to TRUE; the vector is recycled if
#' needed.}
#'
#' \item{move.mu}{a logical indicating whether the mutation rates
#' should be estimated ('moved' in the MCMC), or not, all defaulting to TRUE.}
#'
#' \item{move.pi}{a logical indicating whether the reporting probability
#' should be estimated ('moved' in the MCMC), or not, all defaulting to TRUE.}
#' 
#' \item{move.eps}{a logical indicating whether the contact reporting probability
#' should be estimated ('moved' in the MCMC), or not, all defaulting to TRUE.}
#'
#' \item{move.kappa}{a logical indicating whether the number of generations between two successive
#' cases should be estimated ('moved' in the MCMC), or not, all defaulting to TRUE.}
#'
#' \item{move.pi}{a logical indicating whether the reporting probability
#' should be estimated ('moved' in the MCMC), or not, all defaulting to TRUE.}
#'
#' \item{n.iter}{the number of iterations of the MCMC}
#'
#' \item{sample.every}{the frequency at which MCMC samples are retained for the output}
#'
#' \item{sd.mu}{the standard deviation for the Normal proposal for the mutation rates}
#'
#' \item{sd.pi}{the standard deviation for the Normal proposal for the reporting probability}
#' 
#'  \item{sd.eps}{the standard deviation for the Normal proposal for the contact reporting probability}
#'
#' \item{prop.alpha.move}{the proportion of ancestries to move at each iteration of the MCMC}
#'
#' \item{prop.t.inf.move}{the proportion of infection dates to move at each iteration of the MCMC}

#' \item{batch.size}{the size of the batch of random number pre-generated}
#'
#' \item{paranoid}{a logical indicating if the paranoid mode should be used; this mode is used for
#' performing additional tests during outbreaker; it makes computations substantially slower and is
#' mostly used for debugging purposes.}
#'
#' \item{min.date}{earliest infection date possible, expressed as days since the first sampling}
#'
#' \item{max.kappa}{an integer indicating the largest number of generations between any two linked cases; defaults to 5}
#'
#' \item{prior.mu}{a numeric value indicating the rate of the exponential prior for the mutation rate 'mu'}
#'
#' \item{prior.pi}{a numeric vector of length 2 indicating the first and second parameter of the
#' beta prior for the reporting probability 'pi'}
#' 
#' \item{prior.eps}{a numeric vector of length 2 indicating the upper and lower bounds of the
#' uniform prior for the contact reporting probability 'eps'}
#'
#' }
#'
#' @param ... settings to be passed to outbreaker
#'
#' @param data an optional list of data items as returned by \code{outbreaker.data}; if provided,
#' this allows for further checks of the outbreaker setings.
#'
#' @param config a previous set of settings, as returned by \code{outbreaker.config}
#'
#' @seealso outbreaker.data to check and process data for outbreaker
#'
#' @author Thibaut Jombart (\email{t.jombart@@imperial.ac.uk})
#'
#' @export
#'
#' @examples
#' ## see default settings
#' outbreaker.config()
#'
#' ## change defaults
#' outbreaker.config(move.alpha=FALSE, n.iter=2e5, sample.every=1000)
#'
#'
#' @importFrom utils modifyList
#'
outbreaker.config <- function(..., data=NULL, config=NULL) {
    ## PROCESS ... ONLY IF NO CONFIG IS PASSED
    if (is.null(config)) {
        config <- list(...)
        if (length(config) == 1L && inherits(config, "outbreaker.config")) {
            config <- config[[1]]
        }
    }
    if (inherits(config, "outbreaker.config")) {
        return(config)
    }

    ## SET DEFAULTS ##
    defaults <- list(init.tree=c("seqTrack","star","random"),
                     init.mu=1e-4,
                     init.t.inf=NULL,
                     init.alpha=NULL,
                     init.kappa=1,
                     init.pi=0.9,
                     init.eps=0.9,
                     move.alpha=TRUE, move.swap.cases=TRUE, move.t.inf=TRUE,
                     move.mu=TRUE, move.kappa=TRUE, move.pi=TRUE, move.eps=TRUE,
                     n.iter=1e4, sample.every=200, sd.mu=0.0001, sd.pi=0.01, sd.eps=0.01,
                     prop.alpha.move=1/4,
                     prop.t.inf.move=0.2,
                     batch.size=1e6,
                     paranoid=FALSE,
                     min.date=-10,
                     max.kappa=5,
                     find.import=TRUE,
                     outlier.threshold=5,
                     n.iter.import=5000,
                     sample.every.import=100,
                     prior.mu=1000,
                     prior.pi=c(10,1),
                     prior.eps=c(0.5,1))

    ## MODIFY CONFIG WITH ARGUMENTS ##
    config <- modify.defaults(defaults, config)

    ## CHECK CONFIG ##
    ## check init.tree
    if (is.character(config$init.tree)) {
        config$init.tree <- match.arg(config$init.tree, c("seqTrack","star","random"))
    }
    if (is.numeric(config$init.tree)) {
        config$init.alpha <- config$init.tree
    }

    ## check / process init.t.inf
    if (!is.null(config$init.t.inf)) {
        if (inherits(config$init.t.inf, "Date")) {
            config$init.t.inf <- config$init.t.inf-min(config$init.t.inf)
        }
        if (inherits(config$init.t.inf, "POSIXct")) {
            config$init.t.inf <- difftime(config$init.t.inf, min(config$init.t.inf), units="days")
        }
        config$init.t.inf <- as.integer(round(config$init.t.inf))
    }

    ## check init.mu
    if (!is.numeric(config$init.mu)) {
        stop("init.mu is not a numeric value")
    }
    if (config$init.mu < 0) {
        stop("init.mu is negative")
    }
    if (config$init.mu > 1) {
        stop("init.mu is greater than 1")
    }
    if (!is.finite(config$init.mu)) {
        stop("init.mu is infinite or NA")
    }

    ## check init.kappa
    if (!is.numeric(config$init.kappa)) {
        stop("init.kappa is not a numeric value")
    }
    config$init.kappa <- as.integer(round(config$init.kappa))
    if (any(config$init.kappa < 1)) {
        stop("init.kappa has values smaller than 1")
    }
    if (any(config$init.kappa > config$max.kappa)) {
        config$init.kappa[config$init.kappa > config$max.kappa] <- config$max.kappa
        warning("values of init.kappa greater than max.kappa have been set to max.kappa")
    }
    if (!all(is.finite(config$init.kappa))) {
        stop("init.kappa has values which are infinite or NA")
    }


    ## check init.pi
    if (!is.numeric(config$init.pi)) {
        stop("init.pi is not a numeric value")
    }
    if (config$init.pi < 0) {
        stop("init.pi is negative")
    }
    if (config$init.pi > 1) {
        stop("init.pi is greater than 1")
    }
    if (!is.finite(config$init.pi)) {
        stop("init.pi is infinite or NA")
    }

    ## check init.pi
    if (!is.numeric(config$init.eps)) {
      stop("init.eps is not a numeric value")
    }
    if (config$init.eps < 0) {
      stop("init.eps is negative")
    }
    if (config$init.eps > 1) {
      stop("init.eps is greater than 1")
    }
    if (!is.finite(config$init.eps)) {
      stop("init.eps is infinite or NA")
    }
    
    ## check move.alpha
    if (!is.logical(config$move.alpha)) {
        stop("move.alpha is not a logical")
    }
    if (is.na(config$move.alpha)) {
        stop("move.alpha is NA")
    }

    ## check move.swap.cases
    if (!is.logical(config$move.swap.cases)) {
        stop("move.swap.cases is not a logical")
    }
    if (is.na(config$move.swap.cases)) {
        stop("move.swap.cases is NA")
    }

    ## check move.t.inf
    if (!is.logical(config$move.t.inf)) {
        stop("move.t.inf is not a logical")
    }
    if (is.na(config$move.t.inf)) {
        stop("move.t.inf is NA")
    }

    ## check move.mu
    if (!is.logical(config$move.mu)) {
        stop("move.mu is not a logical")
    }
    if (is.na(config$move.mu)) {
        stop("move.mu is NA")
    }

    ## check move.kappa
    if (!is.logical(config$move.kappa)) {
        stop("move.kappa is not a logical")
    }
    if (is.na(config$move.kappa)) {
        stop("move.kappa is NA")
    }

    ## check move.pi
    if (!is.logical(config$move.pi)) {
        stop("move.pi is not a logical")
    }
    if (is.na(config$move.pi)) {
        stop("move.pi is NA")
    }

    ## check move.pi
    if (!is.logical(config$move.eps)) {
      stop("move.eps is not a logical")
    }
    if (is.na(config$move.eps)) {
      stop("move.eps is NA")
    }
    
    ## check n.iter
    if (!is.numeric(config$n.iter)) {
        stop("n.iter is not a numeric value")
    }
    if (config$n.iter < 2 ) {
        stop("n.iter is smaller than 2")
    }
    if (!is.finite(config$n.iter)) {
        stop("n.iter is infinite or NA")
    }


    ## check sample.every
    if (!is.numeric(config$sample.every)) {
        stop("sample.every is not a numeric value")
    }
    if (config$sample.every < 1 ) {
        stop("sample.every is smaller than 1")
    }
    if (!is.finite(config$sample.every)) {
        stop("sample.every is infinite or NA")
    }

    ## check sd.mu
    if (!is.numeric(config$sd.mu)) {
        stop("sd.mu is not a numeric value")
    }
    if (config$sd.mu < 1e-10) {
        stop("sd.mu is close to zero or negative")
    }
    if (!is.finite(config$sd.mu)) {
        stop("sd.mu is infinite or NA")
    }

    ## check sd.pi
    if (!is.numeric(config$sd.pi)) {
        stop("sd.pi is not a numeric value")
    }
    if (config$sd.pi < 1e-10) {
        stop("sd.pi is close to zero or negative")
    }
    if (!is.finite(config$sd.pi)) {
        stop("sd.pi is infinite or NA")
    }

    ## check sd.pi
    if (!is.numeric(config$sd.eps)) {
      stop("sd.eps is not a numeric value")
    }
    if (config$sd.eps < 1e-10) {
      stop("sd.eps is close to zero or negative")
    }
    if (!is.finite(config$sd.eps)) {
      stop("sd.eps is infinite or NA")
    }
    
    ## check prop.alpha.move
    if (!is.numeric(config$prop.alpha.move)) {
        stop("prop.alpha.move is not a numeric value")
    }
    if (config$prop.alpha.move < 0 ) {
        stop("prop.alpha.move is negative")
    }
    if (config$prop.alpha.move > 1 ) {
        stop("prop.alpha.move is greater than one")
    }
    if (!is.finite(config$prop.alpha.move)) {
        stop("prop.alpha.move is infinite or NA")
    }

    ## check prop.t.inf.move
    if (!is.numeric(config$prop.t.inf.move)) {
        stop("prop.t.inf.move is not a numeric value")
    }
    if (config$prop.t.inf.move < 0 ) {
        stop("prop.t.inf.move is negative")
    }
    if (config$prop.t.inf.move > 1 ) {
        stop("prop.t.inf.move is greater than one")
    }
    if (!is.finite(config$prop.t.inf.move)) {
        stop("prop.t.inf.move is infinite or NA")
    }

    ## check batch.size
    if (!is.numeric(config$batch.size)) {
        stop("batch.size is not numeric")
    }
    if (config$batch.size<1) {
        stop("batch.size is less than 1")
    }
    if (!is.finite(config$batch.size)) {
        stop("batch.size is infinite or NA")
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

    ## check min.date
    if (!is.numeric(config$min.date)) {
        stop("min.date is not numeric")
    }
    if (config$min.date >= 0) {
        stop("min.date is greater or equal to 0")
    }
    if (!is.finite(config$min.date)) {
        stop("min.date is infinite or NA")
    }

    ## check find.import
    if (!is.logical(config$find.import)) {
        stop("find.import is not logical")
    }
    if (length(config$find.import) != 1L) {
        stop("find.import should be a single logical value")
    }
    if (is.na(config$find.import)) {
        stop("find.import is NA")
    }

    ## check outlier.threshold
    if (!is.numeric(config$outlier.threshold)) {
        stop("outlier.threshold is not a numeric value")
    }
    if (any(config$outlier.threshold < 1)) {
        stop("outlier.threshold has values smaller than 1")
    }
    if (!is.finite(config$outlier.threshold)) {
        stop("outlier.threshold is infinite or NA")
    }

    ## check n.iter.import
    if (!is.numeric(config$n.iter.import)) {
        stop("n.iter.import is not a numeric value")
    }
    if (config$n.iter.import < 1000) {
        stop("n.iter is smaller than 1000")
    }
    if (!is.finite(config$n.iter.import)) {
        stop("n.iter.import is infinite or NA")
    }

    ## check sample.every.import
    if (!is.numeric(config$sample.every.import)) {
        stop("sample.every.import is not a numeric value")
    }
    if (config$sample.every.import < 1) {
        stop("sample.every.import is smaller than 1")
    }
    if (!is.finite(config$sample.every.import)) {
        stop("sample.every.import is infinite or NA")
    }

    ## check prior value for mu
    if (!is.numeric(config$prior.mu)) {
        stop("prior.mu is not a numeric value")
    }
    if (config$prior.mu < 0) {
        stop("prior.mu is negative (it should be a rate)")
    }
    if (!is.finite(config$prior.mu)) {
        stop("prior.mu is infinite or NA")
    }

    ## check prior value for pi
    if (!all(is.numeric(config$prior.pi))) {
        stop("prior.pi has non-numeric values")
    }
    if (any(config$prior.pi < 0)) {
        stop("prior.pi has negative values")
    }
    if (length(config$prior.pi)!=2L) {
        stop("prior.pi should be a vector of length 2")
    }
   if (!all(is.finite(config$prior.pi))) {
        stop("prior.pi is has values which are infinite or NA")
    }
    
    ## check prior value for eps
    if (!all(is.numeric(config$prior.eps))) {
      stop("prior.eps has non-numeric values")
    }
    if (any(config$prior.eps < 0)) {
      stop("prior.eps has negative values")
    }
    if (length(config$prior.eps)!=2L) {
      stop("prior.eps should be a vector of length 2")
    }
    if (!all(is.finite(config$prior.eps))) {
      stop("prior.eps is has values which are infinite or NA")
    }
    
    
    ## CHECKS POSSIBLE IF DATA IS PROVIDED ##
    if (!is.null(data)) {
      ## check initial tree
      if (is.character(config$init.tree)) {
        if (config$init.tree=="seqTrack" && is.null(data$dna)) {
          message("Can't use seqTrack initialization with missing DNA sequences; using a star-like tree")
          config$init.tree <- "star"
        }
        
        ## seqTrack init
        if (config$init.tree=="seqTrack") {
                D.temp <- data$D
                D.temp[!data$can.be.ances] <- 1e30
                config$init.alpha <- apply(D.temp,2,which.min)
                config$init.alpha[data$dates==min(data$dates)] <- NA
                config$init.alpha <- as.integer(config$init.alpha)
            } else if (config$init.tree=="star") {
                config$init.alpha <- rep(which.min(data$dates), length(data$dates))
                config$init.alpha[data$dates==min(data$dates)] <- NA
            } else if (config$init.tree=="random") {
                config$init.alpha <- ralpha(data$dates)
            }
        } else { ## if ancestries are provided
            if (length(config$init.alpha) != data$N) {
                stop("inconvenient length for init.alpha")
            }
            unknownAnces <- config$init.alpha<1 | config$init.alpha>data$N
            if (any(stats::na.omit(unknownAnces))) {
                warning("some initial ancestries refer to unknown cases (idx<1 or >N)")
                config$init.alpha[unknownAnces] <- NA
            }
        }

        ## check initial t.inf
        if (!is.null(config$init.t.inf)) {
            if (any(config$init.t.inf >= data$dates)) {
                stop("Initial dates of infection come after sampling dates / dates of onset.")
            }
        }

        ## recycle move.alpha
        config$move.alpha <- rep(config$move.alpha, length.out=data$N)

        ## recycle move.t.inf
        config$move.t.inf <- rep(config$move.t.inf, length.out=data$N)

        ## recycle move.kappa
        config$move.kappa <- rep(config$move.kappa, length.out=data$N)

        ## recycle init.kappa
        config$init.kappa <- rep(config$init.kappa, length.out=data$N)
        config$init.kappa[is.na(config$init.alpha)] <- NA

        ## disable moves for imported cases
        config$move.alpha[is.na(config$init.alpha)] <- FALSE
        config$move.kappa[is.na(config$init.alpha)] <- FALSE
    }

    ## output is a list of checked settings with a dedicated class (for printing)
    class(config) <- c("outbreaker.config", "list")
    return(config)
}

