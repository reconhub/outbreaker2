#' Process input data for outbreaker
#'
#' This function performs various checks on input data given to outbreaker.  It
#' takes a list of named items as input, performs various checks, set defaults
#' where arguments are missing, and return a correct list of data input. If no
#' input is given, it returns the default settings.
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
#' \item{ctd}{the contact tracing data provided as a matrix or dataframe of two
#' columns, indicating a reported contact between the two individuals whose ids
#' are provided in a given row of the data.}
#'
#' \item{w_dens}{a vector of numeric values indicating the generation time
#' distribution, reflecting the infectious potential of a case t = 1, 2, ...
#' time steps after infection. By convention, it is assumed that
#' newly infected patients cannot see new infections on the same time step. If not
#' standardized, this distribution is rescaled to sum to 1.}
#'
#' \item{f_dens}{similar to \code{w_dens}, except that this is the distribution
#' of the colonization time, i_e. time interval during which the pathogen can
#' be sampled from the patient.}
#'
#'}
#'
#' @param ... a list of data items to be processed (see description)
#'
#' @param data optionally, an existing list of data item as returned by \code{outbreaker_data}.
#'
#' @author Thibaut Jombart (\email{thibautjombart@@gmail.com})
#'
#' @export
#'
#' @examples
#'
#' x <- fake_outbreak
#' outbreaker_data(dates = x$sample, dna = x$dna, w_dens = x$w)
#'
outbreaker_data <- function(..., data = list(...)) {

    ## SET DEFAULTS ##
    defaults <- list(dates = NULL, w_dens = NULL, f_dens = NULL,
                     dna = NULL, ctd = NULL, N = 0L, L = 0L, D = NULL,
                     max_range = NA, can_be_ances = NULL,
                     log_w_dens = NULL, log_f_dens = NULL,
                     contacts = NULL, C_combn = NULL, C_nrow = NULL)

    ## MODIFY DATA WITH ARGUMENTS ##
    data <- modify_defaults(defaults, data, FALSE)


    ## CHECK DATA ##
    ## CHECK DATES
    if (!is.null(data$dates)) {
        if (inherits(data$dates, "Date")) {
            data$dates <- data$dates-min(data$dates)
        }
        if (inherits(data$dates, "POSIXct")) {
            data$dates <- difftime(data$dates, min(data$dates), units="days")
        }
        data$dates <- as.integer(round(data$dates))
        data$N <- length(data$dates)
        data$max_range <- diff(range(data$dates))
        ## get temporal ordering constraint:
        ## canBeAnces[i,j] is 'i' can be ancestor of 'j'
        data$can_be_ances <- outer(data$dates,
                                   data$dates,
                                   FUN="<") # strict < is needed as we impose w(0)=0
        diag(data$can_be_ances) <- FALSE
    }

    ## CHECK W_DENS
    if (!is.null(data$w_dens)) {
        if (any(data$w_dens<0)) {
            stop("w_dens has negative entries (these should be probabilities!)")
        }

        if (any(!is.finite(data$w_dens))) {
            stop("non-finite values detected in w_dens")
        }


        ## Remove trailing zeroes to prevent starting with -Inf temporal loglike
        if(data$w_dens[length(data$w_dens)] < 1e-15) {
            final_index <- max(which(data$w_dens > 1e-15))
            data$w_dens <- data$w_dens[1:final_index]
            warning("Removed trailing zeroes found in w_dens")
        }

        ## add an exponential tail summing to 1e-4 to 'w'
        ## to cover the span of the outbreak
        ## (avoids starting with -Inf temporal loglike)
        if (length(data$w_dens) < data$max_range) {
            length_to_add <- (data$max_range-length(data$w_dens)) + 10 # +10 to be on the safe side
            val_to_add <- stats::dexp(seq_len(length_to_add), 1)
            val_to_add <- 1e-4*(val_to_add/sum(val_to_add))
            data$w_dens <- c(data$w_dens, val_to_add)
        }

        ## standardize the mass function
        data$w_dens <- data$w_dens / sum(data$w_dens)
        data$log_w_dens <- matrix(log(data$w_dens), nrow = 1)
    }

    ## CHECK F_DENS
    if (!is.null(data$w_dens) && is.null(data$f_dens)) {
        data$f_dens <- data$w_dens
    }
    if (!is.null(data$f_dens)) {
        if (any(data$f_dens<0)) {
            stop("f_dens has negative entries (these should be probabilities!)")
        }

        if (any(!is.finite(data$f_dens))) {
            stop("non-finite values detected in f_dens")
        }

        data$f_dens <- data$f_dens / sum(data$f_dens)
        data$log_f_dens <- log(data$f_dens)
    }

    
    ## CHECK DNA

    if (!is.null(data$dna)) {
        if (!inherits(data$dna, "DNAbin")) stop("dna is not a DNAbin object.")
        if (!is.matrix(data$dna)) data$dna <- as.matrix(data$dna)

        ## get matrix of distances
        
        data$L <- ncol(data$dna) #  (genome length)
        data$D <- as.matrix(ape::dist.dna(data$dna, model="N")) # distance matrix
        storage.mode(data$D) <- "integer" # essential for C/C++ interface
        
        ## get matching between sequences and cases
        
        if (is.null(rownames(data$dna))) {
            if (nrow(data$dna) != data$N) {
                msg <- sprintf(paste("numbers of sequences and cases differ (%d vs %d):",
                                     "please label sequences"),
                               nrow(data$dna), data$N)
                stop(msg)
            }
                
            rownames(data$dna) <- rownames(data$D) <- colnames(data$D) <- seq_len(data$N)
        }

        data$id_in_dna <- match(as.character(seq_len(data$N)), rownames(data$dna))
        
    } else {
        data$L <- 0L
        data$D <- matrix(integer(0), ncol = 0, nrow = 0)
        data$id_in_dna <- rep(NA_integer_, data$N)
    }
    data$has_dna <- !is.na(data$id_in_dna)
    

    ## CHECK CTD
    
    if (!is.null(data$ctd)) {
        if (!inherits(data$ctd, c("matrix", "data.frame"))) {
            stop("ctd is not a matrix or data.frame")
        }
        if (!is.matrix(data$ctd)) data$ctd <- as.matrix(data$ctd)
        not.found <- data$ctd[any(!data$ctd %in% 1:data$N)]
        if (length(not.found) != 0) {
            not.found <- sort(unique(not.found))
            stop(paste("Individual(s)", paste(not.found, collapse = ", "),
                       "are unknown cases (idx < 1 or > N")
                 )
        }
        contacts <- matrix(0, data$N, data$N)
        for(i in seq_len(nrow(data$ctd))) {
            pair <- data$ctd[i,]
            contacts[pair[[1]], pair[[2]]] <- contacts[pair[[2]], pair[[1]]] <- 1
        }
        data$contacts <- contacts
        data$C_combn <- data$N*(data$N - 1)/2
        data$C_nrow <- nrow(data$ctd)
    } else {
        data$contacts <- matrix(integer(0), ncol = 0, nrow = 0)
    }

    ## output is a list of checked data
    return(data)

}

