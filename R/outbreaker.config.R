#' Set and check parameter settings of outbreaker
#'
#' This function defines settings for outbreaker.
#' It takes a list of named items as input, performs various checks, set defaults where arguments are missing, and return a correct list of settings. If no input is given, it returns the default settings.
#'
#' Acceptables arguments for ... are:
#' \describe{
#' \item{init.tree}{the tree used to initialize the MCMC. Can be either a
#' character string indicating how this tree should be computed, or a vector of
#' integers corresponding to the tree itself, where the i-th value corresponds
#' to the index of the ancestor of 'i' (i.e., \code{init.tree[i]} is the
#' ancestor of case \code{i}). Accepted character strings are "seqTrack" (uses
#' seqTrack output as initialize tree), "random" (ancestor randomly selected
#' from preceding cases), and "star" (all cases coalesce to the first case).
#' Note that for SeqTrack, all cases should have been sequenced.}
#'
#' \item{n.iter}{an integer indicating the number of iterations in the MCMC,
#' including the burnin period}
#'
#' \item{init.mu}{initial values for the mutation rates}
#'
#' \item{move.ances}{a vector of logicals indicating, for each case, if the ancestry should be estimated ('moved' in the MCMC), or not, defaulting to TRUE; the vector is recycled if needed.}
#'
#' \item{move.t.inf}{a vector of logicals indicating, for each case, if the dates of infection should be estimated ('moved' in the MCMC), or not, defaulting to TRUE; the vector is recycled if needed.}
#'
#' \item{move.mu}{a logical indicating whether the mutation rates
#' should be estimated ('moved' in the MCMC), or not, all defaulting to TRUE.}
#'
#' \item{n.iter}{the number of iterations of the MCMC}
#'
#' \item{sample.every}{the frequency at which MCMC samples are retained for the output}
#'
#' \item{sd.mu}{the standard deviation for the Normal proposal for the mutation rates}
#'
#' \item{prop.ances.move}{the proportion of ancestries to move at each iteration of the MCMC}
#' }
#'
#' @param ... settings to be passed to outbreaker
#'
#' @param data an optional list of data items as returned by \code{outbreaker.data}; if provided, this allows for further checks of the outbreaker setings.
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
#' outbreaker.config(move.ances=FALSE, n.iter=2e5, sample.every=1000)
#'
#'
#' @importFrom utils modifyList
#'
outbreaker.config <- function(..., data=NULL, config=NULL){
    ## PROCESS ... ONLY IF NO CONFIG IS PASSED
    if(is.null(config)){
        config <- list(...)
    }

    ## SET DEFAULTS ##
    defaults <- list(init.tree=c("seqTrack","star","random"),
                     init.mu=1e-4,
                     init.t.inf=NULL,
                     move.ances=TRUE, move.t.inf=TRUE, move.mu=TRUE,
                     n.iter=100, sample.every=10, sd.mu=0.0001,
                     prop.ances.move=1/4,
                     batch.size=1e6)

    ## MODIFY CONFIG WITH ARGUMENTS ##
    config <- modify.defaults(defaults, config)

    ## CHECK CONFIG ##
    ## check init.tree
    if(is.character(config$init.tree)){
        config$init.tree <- match.arg(config$init.tree, c("seqTrack","star","random"))
    }

    ## check / process init.t.inf
    if(!is.null(config$init.t.inf)){
        if(inherits(config$init.t.inf, "Date")) config$init.t.inf <- config$init.t.inf-min(config$init.t.inf)
        if(inherits(config$init.t.inf, "POSIXct")) config$init.t.inf <- difftime(config$init.t.inf, min(config$init.t.inf), units="days")
        config$init.t.inf <- as.integer(round(config$init.t.inf))
    }

    ## check init.mu
    if(!is.numeric(config$init.mu)) stop("init.mu is not a numeric value")
    if(config$init.mu < 0) stop("init.mu is negative")
    if(config$init.mu > 1) stop("init.mu is greater than 1")

    ## check move.ances
    if(!is.logical(config$move.ances)) stop("move.ances is not a logical")

    ## check move.t.inf
    if(!is.logical(config$move.t.inf)) stop("move.t.inf is not a logical")

    ## check move.mu
    if(!is.logical(config$move.mu)) stop("move.mu is not a logical")

    ## check n.iter
    if(!is.numeric(config$n.iter)) stop("n.iter is not a numeric value")
    if(config$n.iter < 2 ) stop("n.iter is smaller than 2")

    ## check sample.every
    if(!is.numeric(config$sample.every)) stop("sample.every is not a numeric value")
    if(config$sample.every < 1 ) stop("sample.every is smaller than 1")

    ## check sd.mu
    if(!is.numeric(config$sd.mu)) stop("sd.mu is not a numeric value")
    if(config$sd.mu < 0) stop("sd.mu is negative")

    ## check prop.ances.move
    if(!is.numeric(config$prop.ances.move)) stop("prop.ances.move is not a numeric value")
    if(config$prop.ances.move < 0 ) stop("prop.ances.move is negative")
    if(config$prop.ances.move > 1 ) stop("prop.ances.move is greater than one")

    ## check batch.size
    if(!is.numeric(config$batch.size)) stop("batch.size is not numeric")
    if(config$batch.size<1) stop("batch.size is less than 1")


    ## CHECKS POSSIBLE IF DATA IS PROVIDED ##
    if(!is.null(data)){
        ## check initial tree
        if(is.character(config$init.tree)){
            if(config$init.tree=="seqTrack" && is.null(data$dna)){
                warning("Can't use seqTrack initialization with missing DNA sequences; using a star-like tree")
                config$init.tree <- "star"
            }

            ## seqTrack init
            if(config$init.tree=="seqTrack"){
                D.temp <- data$D
                D.temp[!data$CAN.BE.ANCES] <- 1e30
                config$ances <- apply(D.temp,2,which.min)
                config$ances[data$dates==min(data$dates)] <- NA
                config$ances <- as.integer(config$ances)
            } else if(config$init.tree=="star"){
                config$ances <- rep(which.min(data$dates), length(data$dates))
                config$ances[data$dates==min(data$dates)] <- NA
            } else if(config$init.tree=="random"){
                config$ances <- rances(data$dates)
            }
        } else { ## if ancestries are provided
            if(length(config$init.tree) != data$N) stop("inconvenient length for init.tree")
            unknownAnces <- config$init.tree<1 | config$init.tree>data$N
            if(any(na.omit(unknownAnces))){
                warning("some initial ancestries refer to unknown cases (idx<1 or >N)")
                config$init.tree[unknownAnces] <- NA
            }
        }

        ## check initial t.inf
        if(!is.null(config$init.t.inf)){
            if(any(config$init.t.inf >= data$dates)) stop("Initial dates of infection come after sampling dates / dates of onset.")
        }

        ## recycle move.ances
        config$move.ances <- rep(config$move.ances, lenth=data$N)

        ## recycle move.t.inf
        config$move.t.inf <- rep(config$move.t.inf, lenth=data$N)
    } else {
        ## set ances to NULL
        config$ances <- NULL
    }

    ## RETURN CONFIG ##
    return(config)
} # end outbreaker.config





## NON EXPORTED FUNCTIONS ##

## append arguments to defaults ##
modify.defaults <- function(defaults, x){
    extra <- setdiff(names(x), names(defaults))
    if (length(extra) > 0L){
        stop("Additional invalid options: ", paste(extra, collapse=", "))
    }
    modifyList(defaults, x)
} # end modify.defaults
