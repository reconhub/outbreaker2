#' Checks input data for outbreaker
#'
#' This function performs various checks on input data given to outbreaker.
#'
#' @param dna the DNA sequences in \code{DNAbin} format (see
#' \code{\link[ape]{read.dna}} in the ape package); this can be imported from a
#' fasta file (extension .fa, .fas, or .fasta) using \code{adegenet}'s function
#' \link[adegenet]{fasta2DNAbin}.
#'
#' @param dates a vector indicating the collection dates, provided either as
#' integer numbers or in a usual date format such as \code{Date} or
#' \code{POSIXct} format. By convention, zero will indicate the oldest date.
#'
#' @param w.dens a vector of numeric values indicating the generation time
#' distribution, reflecting the infectious potential of a case t=0, 1, 2, ...
#' time steps after infection. By convention, w.dens[1]=0, meaning that an
#' newly infected patient cannot be instantaneously infectious. If not
#' standardized, this distribution is rescaled to sum to 1.
#'
#' @param f.dens similar to \code{w.dens}, except that this is the distribution
#' of the colonization time, i.e. time interval during which the pathogen can
#' be sampled from the patient.
#'
#' @author Thibaut Jombart (\email{t.jombart@@imperial.ac.uk})
#'
#' @export
#'
#' @examples
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
    ## get temporal ordering constraint:
    ## canBeAnces[i,j] is 'i' can be ancestor of 'j'
    CAN.BE.ANCES <- outer(dates,dates,FUN="<") # strict < is needed as we impose w(0)=0
    diag(CAN.BE.ANCES) <- FALSE


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
                N=N, L=L, D=D, MAX.RANGE=MAX.RANGE, CAN.BE.ANCES=CAN.BE.ANCES,
                log.w.dens=log.w.dens, log.f.dens=log.f.dens))

} # end outbreaker.data





