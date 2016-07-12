#' Process input data for outbreaker
#'
#' This function performs various checks on input data given to outbreaker.
#' It takes a list of named items as input, performs various checks, set defaults where arguments are missing, and return a correct list of data input. If no input is given, it returns the default settings.
#'
#' Acceptables arguments for ... are:
#' \describe{
#' \item{dates}{dates a vector indicating the collection dates, provided either as
#' integer numbers or in a usual date format such as \code{Date} or
#' \code{POSIXct} format. By convention, zero will indicate the oldest date.}
#'
#' \item{dna}{the DNA sequences in \code{DNAbin} format (see
#' \code{\link[ape]{read.dna}} in the ape package); this can be imported from a
#' fasta file (extension .fa, .fas, or .fasta) using \code{adegenet}'s function
#' \link[adegenet]{fasta2DNAbin}.}
#' 
#' \item{CTD}{the contact tracing data provided as a matrix or dataframe of two columns,
#' containing the pairs of id's with observed contact}
#'
#' \item{w.dens}{a vector of numeric values indicating the generation time
#' distribution, reflecting the infectious potential of a case t=1, 2, ...
#' time steps after infection. By convention, it is assumed that
#' newly infected patients cannot see new infections on the same time step. If not
#' standardized, this distribution is rescaled to sum to 1.}
#'
#' \item{f.dens}{similar to \code{w.dens}, except that this is the distribution
#' of the colonization time, i.e. time interval during which the pathogen can
#' be sampled from the patient.}
#'
#'}
#'
#' @param ... a list of data items to be processed (see description)
#'
#' @param data optionally, an existing list of data item as returned by \code{outbreaker.data}.
#'
#' @author Thibaut Jombart (\email{t.jombart@@imperial.ac.uk})
#'
#' @export
#'
#' @examples
#'
#' x <- fake.outbreak
#' outbreaker.data(dates=x$collecDates, dna=x$dat$dna, w.dens=x$w)
#'
outbreaker.data <- function(..., data=list(...)) {

    ## SET DEFAULTS ##
    defaults <- list(dates=NULL, w.dens=NULL, f.dens=NULL, dna=NULL, CTD=NULL, contact=NULL,
                     N=0L, L=0L, D=NULL, max.range=NA, can.be.ances=NULL,
                     log.w.dens=NULL, log.f.dens=NULL)

    ## MODIFY DATA WITH ARGUMENTS ##
    data <- modify.defaults(defaults, data)


    ## CHECK DATA ##
    ## CHECK DATES
    if (!is.null(data$dates)) {
        if (inherits(data$dates, "Date")) data$dates <- data$dates-min(data$dates)
        if (inherits(data$dates, "POSIXct")) data$dates <- difftime(data$dates, min(data$dates), units="days")
        data$dates <- as.integer(round(data$dates))
        data$N <- length(data$dates)
        data$max.range <- diff(range(data$dates))
        ## get temporal ordering constraint:
        ## canBeAnces[i,j] is 'i' can be ancestor of 'j'
        data$can.be.ances <- outer(data$dates,data$dates,FUN="<") # strict < is needed as we impose w(0)=0
        diag(data$can.be.ances) <- FALSE
    }

    ## CHECK W.DENS
    if (!is.null(data$w.dens)) {
        if (any(data$w.dens<0))
        {
            stop("w.dens has negative entries (these should be probabilities!)")
        }
        ## add an exponential tail summing to 1e-4 to 'w'
        ## to cover the span of the outbreak
        ## (avoids starting with -Inf temporal loglike)
        if (length(data$w.dens)<data$max.range) {
            length.to.add <- (data$max.range-length(data$w.dens)) + 10 # +10 to be on the safe side
            val.to.add <- stats::dexp(seq_len(length.to.add), 1)
            val.to.add <- 1e-4*(val.to.add/sum(val.to.add))
            data$w.dens <- c(data$w.dens, val.to.add)
            data$w.dens <- data$w.dens/sum(data$w.dens)
        }
        data$log.w.dens <- matrix(log(data$w.dens), nrow=1)
    }

    ## CHECK F.DENS
    if (!is.null(data$w.dens) && is.null(data$f.dens)) {
        data$f.dens <- data$w.dens
    }
    if (!is.null(data$f.dens)) {
        if (any(data$f.dens<0))
        {
            stop("f.dens has negative entries (these should be probabilities!)")
        }
        data$log.f.dens <- log(data$f.dens)
    }

    ## CHECK DNA
    if (!is.null(data$dna)) {
        if (!inherits(data$dna, "DNAbin")) stop("dna is not a DNAbin object.")
        if (!is.matrix(data$dna)) data$dna <- as.matrix(data$dna)
        data$L <- ncol(data$dna) #  (genome length)
        data$D <- as.matrix(ape::dist.dna(data$dna, model="N")) # distance matrix
    } else {
        data$L <- 0
        data$D <- matrix(numeric(0), ncol=0, nrow=0)
    }
    
    ## CHECK CONTACT
    if (!is.null(data$CTD)) {
      contact <- matrix(FALSE,data$N,data$N)
      for(cij in seq_len(nrow(data$CTD))) contact[data$CTD[cij,1],data$CTD[cij,2]] <- TRUE
      data$contact <- contact | t(contact)
    } else {
      data$contact <- matrix(numeric(0), ncol=0, nrow=0)
    }

    ## output is a list of checked data
    return(data)

}





