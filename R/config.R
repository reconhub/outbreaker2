#' Checks input data for outbreaker
#'
#' This function performs various checks on input data given to outbreaker.
#'
#' Acceptables arguments are:
#' \describe{
#'
#' \item{dna}{the DNA sequences in \code{DNAbin} format (see
#' \code{\link[ape]{read.dna}} in the ape package); this can be imported from a
#' fasta file (extension .fa, .fas, or .fasta) using \code{adegenet}'s function
#' \link[adegenet]{fasta2DNAbin}.}
#'
#' \item{dates}{a vector indicating the collection dates, provided either as
#' integer numbers or in a usual date format such as \code{Date} or
#' \code{POSIXct} format. By convention, zero will indicate the oldest date.}
#'
#' \item{w.dens}{a vector of numeric values indicating the generation time
#' distribution, reflecting the infectious potential of a case t=0, 1, 2, ...
#' time steps after infection. By convention, w.dens[1]=0, meaning that an
#' newly infected patient cannot be instantaneously infectious. If not
#' standardized, this distribution is rescaled to sum to 1.}
#'
#' \item{f.dens}{similar to \code{w.dens}, except that this is the distribution
#' of the colonization time, i.e. time interval during which the pathogen can
#' be sampled from the patient.}
#'
#' \item{init.tree}{the tree used to initialize the MCMC. Can be either a
#' character string indicating how this tree should be computed, or a vector of
#' integers corresponding to the tree itself, where the i-th value corresponds
#' to the index of the ancestor of 'i' (i.e., \code{init.tree[i]} is the
#' ancestor of case \code{i}). Accepted character strings are "seqTrack" (uses
#' seqTrack output as initialize tree), "random" (ancestor randomly selected
#' from preceding cases), and "star" (all cases coalesce to the first case).
#' Note that for SeqTrack, all cases should have been sequenced.}
#'
#'}
#'
#' @author Thibaut Jombart (\email{t.jombart@@imperial.ac.uk})
#'
#' @export
#'
#' x <- fakeOutbreak
#' outbreaker.data(dates=x$collecDates, dna=x$dat$dna, w.dens=x$w)
#'
outbreaker.data <- function(dates, w.dens, f.dens=w.dens, dna=NULL) {
    ## CHECK DATES ##
    if(inherits(dates, "Date")) dates <- dates-min(dates)
    if(inherits(dates, "POSIXct")) dates <- difftime(dates, min(dates), units="days")
    dates <- as.integer(round(dates))
    N <- length(dates)
    MAX.RANGE <- diff(range(dates))


    ## CHECK W.DENS ##
    if(any(w.dens<0)) {
        stop("w.dens has negative entries (these should be probabilities!)")
    }
    w.dens[1] <- 0
    ## add an exponential tail summing to 1e-4 to 'w'
    ## to cover the span of the outbreak
    ## (avoids starting with -Inf temporal loglike)
    if(length(w.dens)<MAX.RANGE) {
        length.to.add <- (MAX.RANGE-length(w.dens)) + 10 # +10 to be on the safe side
        val.to.add <- dexp(1:length.to.add, 1)
        val.to.add <- 1e-4*(val.to.add/sum(val.to.add))
        w.dens <- c(w.dens, val.to.add)
        w.dens <- w.dens/sum(w.dens)
    }
    log.w.dens <- log(w.dens)


    ## CHECK F.DENS ##
    if(any(f.dens<0)) {
        stop("f.dens has negative entries (these should be probabilities!)")
    }
    f.dens[1] <- 0
        log.f.dens <- log(f.dens)


    ## CHECK DNA ##
    if(!is.null(dna)){
        if(!inherits(dna, "DNAbin")) stop("dna is not a DNAbin object.")
        if(!is.matrix(dna)) dna <- as.matrix(dna)
        L <- ncol(dna) #  (genome length)
        D <- as.matrix(dist.dna(dna, model="N")) # distance matrix
    } else {
        L <- 0
        D <- matrix(numeric(0), ncol=0, nrow=0)
    }

    ## RETURN DATA ##
    return(list(dates=dates, dna=dna, w.dens=w.dens, f.dens=f.dens,
                N=N, L=L, D=D, MAX.RANGE=MAX.RANGE,
                log.w.dens=log.w.dens, log.f.dens=log.f.dens))

} # end outbreaker.data






#' #' Set and check parameter settings of outbreaker
#'
#' This function defines settings for outbreaker.
#' It takes a list of named items as input, performs various checks, set defaults where arguments are missing, and return a correct list of settings. If no input is given, it returns the default settings.
#'
#' Acceptables arguments are:
#' \describe{
#'
#' \item{data}{a list of data items as returned by \code{outbreaker.data}}
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
#'}
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
outbreaker.config <- function(data, ...) {
    config <- list(...)

    ## SET DEFAULTS ##
    defaults <- list(init.tree=c("seqTrack","star","random"),
                     init.mu=1e-4,
                     move.ances=TRUE, move.t.inf=TRUE, move.mu=TRUE,
                     n.iter=100, sample.every=10, sd.mu=0.0001)

    ## MODIFY CONFIG WITH ARGUMENTS ##
    config <- modify.defaults(defaults, config)

    ## CHECK CONFIG ##
    ## check init.tree
    config$init.tree <- match.arg(config$init.tree, c("seqTrack","star","random"))

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

    ## RETURN CONFIG ##
    return(config)
} # end outbreaker.config



#' This function is not exported
#'
modify.defaults <- function(defaults, x) {
    extra <- setdiff(names(x), names(defaults))
    if (length(extra) > 0L) {
        stop("Additional invalid options: ", paste(extra, collapse=", "))
    }
    modifyList(defaults, x)
} # end modify.defaults





