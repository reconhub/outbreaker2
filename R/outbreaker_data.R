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
#' \code{POSIXct} format. By convention, zero will indicate the oldest date. If
#' the vector is named, the vector names will be used for matching cases to
#' contact tracing data and labelled DNA sequences.}
#'
#' \item{dna}{the DNA sequences in \code{DNAbin} format (see
#' \code{\link[ape]{read.dna}} in the ape package); this can be imported from a
#' fasta file (extension .fa, .fas, or .fasta) using \code{adegenet}'s function
#' \link[adegenet]{fasta2DNAbin}.}
#'
#' \item{ctd}{the contact tracing data provided as a matrix/dataframe of two
#' columns, indicating a reported contact between the two individuals whose ids
#' are provided in a given row of the data, or an epicontacts object. In the case
#' of the latter, linelist IDs will be used for matching dates and DNA
#' sequences.}
#'
#' \item{w_dens}{a vector of numeric values indicating the generation time
#' distribution, reflecting the infectious potential of a case t = 1, 2, ...
#' time steps after infection. By convention, it is assumed that
#' newly infected patients cannot see new infections on the same time step. If not
#' standardized, this distribution is rescaled to sum to 1.}
#'
#' \item{f_dens}{similar to \code{w_dens}, except that this is the distribution
#' of the colonization time, i_e. time interval during which the pathogen can
#' be sampled from the patient.}}
#'
#' @param ... a list of data items to be processed (see description).
#'
#' @param data optionally, an existing list of data item as returned by \code{outbreaker_data}.
#'
#' @author Thibaut Jombart (\email{thibautjombart@@gmail.com}).
#'
#' @importFrom magrittr %>%
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
  defaults <- list(dates = NULL,
                   w_dens = NULL,
                   f_dens = NULL,
                   dna = NULL,
                   phydna = NULL,
                   ctd = NULL,
                   ctd_timed = NULL,
                   ids = NULL,
                   w_unobs = NULL,
                   N = 0L,
                   L = 0L,
                   D = NULL,
                   max_range = NA,
                   can_be_ances = NULL,
                   log_w_dens = NULL,
                   log_w_unobs = NULL,
                   log_f_dens = NULL,
                   contacts = NULL,
                   C_combn = NULL,
                   C_nrow = NULL,
                   C_ind = NULL,
                   has_dna = logical(0),
                   has_dna_ind = numeric(),
                   dna_combn = NULL,
                   N_place = NULL,
                   N_place_unobserved = NULL,
                   p_trans = NULL,
                   p_place = NULL,
                   pp_trans = NULL,
                   pp_place = NULL,
                   pp_place_adj = NULL,
                   prop_place_observed = NULL,
                   has_ctd_timed = NULL,
                   p_wrong = 0,
                   id_in_f = integer(0),
                   id_in_dna = integer(0))

  ## MODIFY DATA WITH ARGUMENTS ##
  data <- modify_defaults(defaults, data, FALSE)

  ## Set up case ids
  if(is.null(data$ids)) {
    if(!is.null(names(data$dates))) {
      data$ids <- names(data$dates)
    } else if(!is.null(data$ctd) && inherits(data$ctd, "epicontacts")){
      data$ids <- as.character(data$ctd$linelist$id)
    } else {
      data$ids <- as.character(seq_along(data$dates))
    }
  }

  ## CHECK DATA ##
  ## CHECK DATES
  if (!is.null(data$dates)) {
    min_date <- min(data$dates)
    if (inherits(data$dates, "Date")) {
      data$dates <- data$dates - min_date
    }
    if (inherits(data$dates, "POSIXct")) {
      data$dates <- difftime(data$dates, min_date, units="days")
    }
    data$dates <- setNames(as.integer(round(data$dates)), names(data$dates))
    data$dates <- data$dates - min(data$dates)
    ## store the first date for downstream use (e.g. converting t_inf back)
    if (inherits(min_date, "Date") || inherits(min_date, "POSIXct")) {
      data$min_date <- as.Date(min_date)
    }
    ## assign 'default' name for mapping to f_dens
    if(is.null(names(data$dates))) {
      names(data$dates) <- rep("default", length(data$dates))
    }
    if (inherits(data$dates, "numeric") && any(data$dates %% 1 != 0)) {
      warning("Rounding non-integer dates to nearest integer")
    }
    data$N <- length(data$dates)
    data$max_range <- diff(range(data$dates))
  }

  ## CHECK ID
  if(!is.null(data$ids)) {
    data$ids <- as.character(data$ids)
  } else {
    data$ids <- as.character(seq_len(data$N))
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
    if(data$w_dens[length(data$w_dens)] == 0) {
      final_index <- max(which(data$w_dens > 0))
      data$w_dens <- data$w_dens[1:final_index]
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

  ## CHECK W_UNOBS
  if (!is.null(data$w_unobs)) {
    if (any(data$w_unobs<0)) {
      stop("w_unobs has negative entries (these should be probabilities!)")
    }

    if (any(!is.finite(data$w_unobs))) {
      stop("non-finite values detected in w_unobs")
    }

    ## Remove trailing zeroes to prevent starting with -Inf temporal loglike
    if(data$w_unobs[length(data$w_unobs)] == 0) {
      final_index <- max(which(data$w_unobs > 0))
      data$w_unobs <- data$w_unobs[1:final_index]
      warning("Removed trailing zeroes found in w_unobs")
    }

    ## add an exponential tail summing to 1e-4 to 'w'
    ## to cover the span of the outbreak
    ## (avoids starting with -Inf temporal loglike)
    if (length(data$w_unobs) < data$max_range) {
      length_to_add <- (data$max_range-length(data$w_unobs)) + 10 # +10 to be on the safe side
      val_to_add <- stats::dexp(seq_len(length_to_add), 1)
      val_to_add <- 1e-4*(val_to_add/sum(val_to_add))
      data$w_unobs <- c(data$w_unobs, val_to_add)
    }

    ## standardize the mass function
    data$w_unobs <- data$w_unobs / sum(data$w_unobs)
    data$log_w_unobs <- matrix(log(data$w_unobs), nrow = 1)
  } else {
    data$w_unobs <- data$w_dens
    data$log_w_unobs <- data$log_w_dens
  }

  ## CHECK F_DENS
  if (!is.null(data$w_dens) && is.null(data$f_dens)) {
    data$f_dens <- data$w_dens
  }
  if (!is.null(data$f_dens)) {
    if (!inherits(data$f_dens, c("matrix", "numeric"))) {
      stop("f_dens must be a numeric matrix or vector")
    }
    if (!is.matrix(data$f_dens)) {
      data$f_dens <- matrix(data$f_dens, ncol = 1, dimnames = list(NULL, "default"))
    }
    if (any(data$f_dens<0)) {
      stop("f_dens has negative entries (these should be probabilities!)")
    }
    if (any(!is.finite(data$f_dens))) {
      stop("non-finite values detected in f_dens")
    }
    missing_f <- names(data$dates)[!names(data$dates) %in% colnames(data$f_dens)]
    if (length(missing_f) > 0) {
      stop(sprintf(paste(
        "Named date vectors are used to specify incubation periods for",
        "different groups: either provide a matrix of incubation periods",
        "with column names matching each group, or provide an unnamed date",
        "vector. Incubation periods are missing for the following groups: %s"),
        paste0(unique(missing_f), collapse = ", "))
        )
    }
    data$f_dens <- apply(data$f_dens, 2, function(x) x/sum(x))
    data$log_f_dens <- apply(data$f_dens, 2, log)
    data$id_in_f <- match(names(data$dates), colnames(data$f_dens))
  }

  ## CHECK POTENTIAL ANCESTRIES
  # don't recalculate if can_be_ances is manually provided
  if (is.null(data$can_be_ances) && !is.null(data$dates)) {
    ## get temporal ordering constraint:
    ## canBeAnces[i,j] is 'i' can be ancestor of 'j'
    ## Calculate the serial interval from w_dens and f_dens
    .get_SI <- function(w_dens, f_dens) {
      wf <- stats::convolve(w_dens, rev(f_dens), type = 'open')
      conv <- stats::convolve(rev(f_dens), rev(wf), type = 'open')
      lf <- length(f_dens)
      lw <- length(w_dens)
      return(data.frame(x = (-lf + 2):(lw + lf - 1), d = conv))
    }
    ## Check if difference in sampling dates falls within serial interval
    ## This allows for i to infect j even if it sampled after (SI < 0)
    .can_be_ances <- function(date1, date2, SI) {
      tdiff <- date2 - date1
      out <- sapply(tdiff, function(i) return(i %in% SI$x))
      return(out)
    }
    SI <- .get_SI(data$w_dens, data$f_dens[,1])
    data$can_be_ances <- outer(data$dates,
                               data$dates,
                               FUN=.can_be_ances,
                               SI = SI) # strict < is needed as we impose w(0)=0
    diag(data$can_be_ances) <- FALSE
  }

  ## CHECK DNA

  if (!is.null(data$dna)) {
    if (!inherits(data$dna, "DNAbin")) stop("dna is not a DNAbin object.")
    if (!is.matrix(data$dna)) data$dna <- as.matrix(data$dna)

    ## get matrix of distances
    data$L <- ncol(data$dna) #  (genome length)
    if(is.null(data$D)) {
      data$D <- as.matrix(
        ape::dist.dna(data$dna, model="N", pairwise.deletion = TRUE)
      )
    }
    storage.mode(data$D) <- "integer" # essential for C/C++ interface

    ## get matching between sequences and cases

    if (is.null(rownames(data$dna))) {
      if (nrow(data$dna) != data$N) {
        msg <- sprintf(paste("numbers of sequences and cases differ (%d vs %d):",
                             "please label sequences"),
                       nrow(data$dna), data$N)
        stop(msg)
      }

      ## These need to be indices
      rownames(data$D) <- colnames(data$D) <- seq_len(data$N)

      ## These need to match dates/ctd ids
      rownames(data$dna) <- data$ids
    }
    data$id_in_dna <- match(data$ids, rownames(data$dna))
    if(any(is.na(match(rownames(data$dna), data$ids)))) {
      stop("DNA sequence labels don't match case ids")
    }

  } else {
    data$L <- 0L
    data$D <- matrix(integer(0), ncol = 0, nrow = 0)
    data$id_in_dna <- rep(NA_integer_, data$N)
  }
  data$has_dna <- !is.na(data$id_in_dna)

  # get IDs pairs with sequence data and the distance between them
  if(sum(data$has_dna) > 1) {
    data$has_dna_ind <- which(data$has_dna)
    data$dna_combn <- t(combn(data$has_dna_ind, 2))
    # add distance (D rows follow sequence order, not case index 1:N)
    ri <- data$id_in_dna[data$dna_combn[, 1]]
    rj <- data$id_in_dna[data$dna_combn[, 2]]
    if (anyNA(ri) || anyNA(rj)) {
      stop("Internal error: dna_combn references cases without sequence data")
    }
    data$dna_combn <- cbind(
      data$dna_combn,
      data$D[cbind(ri, rj)]
    )
  } else {
    data$dna_combn <- matrix(0, 0, 0)
  }

  ## check DNA dates if provide
  if (!is.null(data$dna_dates)) {
    if (is.null(data$dna)) {
      stop("DNA dates provided but no sequence data")
    } else if (length(data$dna_dates) != nrow(data$dna)) {
      stop(
        sprintf(
          "Different number of dna sequences and dna dates provided (%i vs %i)",
          nrow(data$dna), length(data$dna_dates))
      )
    }
    if (inherits(data$dna_dates, "Date")) {
      data$dna_dates <- data$dna_dates - min_date
    } else if (inherits(data$dna_dates, "POSIXct")) {
      data$dna_dates <- difftime(data$dna_dates, min_date, units = "days")
    } else if (inherits(data$dna_dates, "numeric")) {
      data$dna_dates <- data$dna_dates - min_date
    }
    data$dna_dates <- as.integer(round(data$dna_dates))
  } else {
    if (!is.null(data$dna))
      data$dna_dates <- data$dates[which(!is.na(data$id_in_dna))]
  }


  ## CHECK CTD
  if (!is.null(data$ctd)) {

    # extract contact info from epicontacts object
    if (inherits(data$ctd, "epicontacts"))
      data$ctd <- data$ctd$contacts[
        intersect(c("from", "to", "type"), names(data$ctd$contacts))
      ]

    # convert tbls to dataframes for subsetting
    if (inherits(data$ctd, "tbl")) data$ctd <- as.data.frame(data$ctd)

    if (!inherits(data$ctd, c("matrix", "data.frame")))
      stop("ctd is not a matrix, data.frame or epicontacts object")

    if (!ncol(data$ctd) %in% c(2, 3))
      stop("ctd must have two or three columns")

    ## Convert to character to prevent factors from interfering
    data$ctd[,1] <- as.character(data$ctd[, 1])
    data$ctd[,2] <- as.character(data$ctd[, 2])

    ## Ensure all cases found in linelist
    not_found <- setdiff(unlist(data$ctd[,1:2]), data$ids)
    if (length(not_found) > 0) {
      stop(paste("Individual(s)", paste(not_found, collapse = ", "),
                 "in contact data are found in case IDs"))
    }

    ## Determine the number of contact types
    if (ncol(data$ctd) == 3) {
      if (!inherits(data$ctd[, 3], "factor")) {
        data$ctd[, 3] <- factor(data$ctd[, 3], levels = unique(data$ctd[, 3]))
      }
    } else {
      data$ctd$type <- factor("foo")
    }

    # generate pairwise binary contact matrices, indexed to data$ids
    # rows = potential infectee, cols = potential infector (see cpp_ll_contact)
    data$ctd_matrix <- lapply(
      levels(data$ctd[, 3]),
      function(contact_type) {
        mat <- matrix(0, nrow = data$N, ncol = data$N)
        sub <- data$ctd[data$ctd[, 3] == contact_type, 1:2, drop = FALSE]
        ridx <- match(sub[, 2], data$ids)
        cidx <- match(sub[, 1], data$ids)
        if (anyNA(ridx) || anyNA(cidx)) {
          stop("Contact IDs not found in data$ids when building ctd_matrix")
        }
        mat[cbind(ridx, cidx)] <- 1
        return(mat)
      }
    )

    ## Count the number of each type of contact - these are ordered by their
    ## factor order, which we will use throughout the analysis
    data$C_nrow <- as.vector(table(data$ctd[,3]))

  } else {
    data$ctd_matrix <- list()
  }

  if (!is.null(data$ctd_timed)) {

    if (!inherits(data$ctd_timed, c("matrix", "data.frame"))) {
      stop("ctd_timed is not a matrix or data.frame")
    }
    if (inherits(data$ctd_timed, "tbl")) data$ctd_timed <- as.data.frame(data$ctd_timed)
    if(!ncol(data$ctd_timed) %in% c(4, 5)) {
      stop(paste0("Timed contact data must have four columns (ID | Place",
                  " | Start date | End date), with an optional fifth column",
                  " specifying the type of contact"))
    }
    if(class(min_date) != class(data$ctd_timed[,3]) |
       class(min_date) != class(data$ctd_timed[,4])) {
      stop("Sampling dates and contact dates are not in the same format.")
    }
    if(inherits(data$ctd_timed[,1], "factor")) {
      data$ctd_timed[,1] <- as.character(data$ctd_timed[,1])
    }
    if(!inherits(data$ctd_timed[,1], c("numeric", "character", "integer"))) {
      stop("IDs in contact data must be numbers or characters")
    }
    if(any(!data$ctd_timed[,1] %in% data$ids)) {
      stop("IDs in timed contact data not found in case IDs.")
    }

    ## Replace IDs with their numeric indices
    data$ctd_timed[,1] <- as.character(data$ctd_timed[,1])

    ## Set place names as characters to prevent accidental indexing with numbers
    data$ctd_timed[,2] <- as.character(data$ctd_timed[,2])

    ## Set min(data$date) as day = 0 and adjust contact dates accordingly
    if (inherits(data$ctd_timed[,3], "Date")) {
      data$ctd_timed[,3] <- data$ctd_timed[,3] - min_date
      data$ctd_timed[,4] <- data$ctd_timed[,4] - min_date
    } else if (inherits(data$ctd_timed[,3], "POSIXct")) {
      data$ctd_timed[,3] <- difftime(data$ctd_timed[,3], min_date, units="days")
      data$ctd_timed[,4] <- difftime(data$ctd_timed[,4], min_date, units="days")
    } else if (inherits(data$ctd_timed[,3], "numeric")) {
      data$ctd_timed[,3] <- data$ctd_timed[,3] - min_date
      data$ctd_timed[,4] <- data$ctd_timed[,4] - min_date
    }
    data$ctd_timed[,3] <- as.integer(round(data$ctd_timed[,3]))
    data$ctd_timed[,4] <- as.integer(round(data$ctd_timed[,4]))

    if(any(data$ctd_timed[,4] - data$ctd_timed[,3] < 0)) {
      stop("End dates must be after start dates.")
    }
    if(ncol(data$ctd_timed) == 5) {
      data$ctd_timed[,5] <- factor(data$ctd_timed[,5],
                                   levels = unique(data$ctd_timed[,5]))
    } else {
      data$ctd_timed$type <- factor('foo')
    }

    contact_types <- levels(data$ctd_timed$type)
    n_types <- length(contact_types)

    ## check p_place
    if(!is.null(data$p_place)) {
      if(length(data$p_place) != n_types) {
        stop("p_place must be provided for each timed contact type")
      }
    } else {
      data$p_place <- vector("list", n_types)
    }

    ## check p_trans
    if(!is.null(data$p_trans)) {
      if(length(data$p_trans) != n_types) {
        stop("p_trans must be provided for each timed contact type")
      }
    } else {
      data$p_trans <- vector("list", n_types)
    }

    ## check prop_place_observed
    if(!is.null(data$prop_place_observed)) {
      if(length(data$prop_place_observed) != n_types) {
        stop("prop_place_observed must be provided for each timed contact type")
      }
      if(any(data$prop_place_observed < 0) | any(data$prop_place_observed > 1)) {
        stop("prop_place_observed must be greater than 0 and smaller than or equal 1")
      }
    } else {
      ## default assumption is that we observed all places
      data$prop_place_observed = rep(1, n_types)
    }

    ## list of timelines for each contact type
    data$ctd_timed_matrix <-
      data$pp_place <-
        data$pp_place_adj <-
          data$pp_trans <-
            data$pp_trans_adj <- list()

    ## indexing to go from dates to column index
    data$C_ind <- -min(data$ctd_timed[,3])

    ## the inferred number of unobserved places
    data$N_place_unobserved = numeric(n_types)

    ## pass places as numbers rather than strings to speed up comparison
    ## these numbers will refer to the indices in the transition matrices
    place_order <- tapply(
      data$ctd_timed[,2],
      data$ctd_timed[,5],
      function(x) unique(x[!is.na(x)])
    )

    ## total range of dates
    t_range <- max(data$ctd_timed[,4]) - min(data$ctd_timed[,3]) + 1

    ## this will calculate the pairwise transition probabilities from a vector
    ## of individual probabilties - we specify p(i to i) = 0
    get_trans <- function(x, loc_vec) {
      if(x[2] != x[1]) loc_vec[x[2]]/(1 - loc_vec[x[1]]) else 0
    }

    ## infers the number of unobserved places, given the number of observed
    ## places and the proportion of total places they represent
    get_n_unobs <- function(n_obs, prop) {
      ceiling(n_obs*((1 - prop)/prop))
    }

    ## insert places into timelines - note this will overwrite previous
    ## locations if multiple locations are provided for the same date - we also
    ## calculate transition probabilities for each contact type
    for(i in 1:n_types) {

      sub <- subset(data$ctd_timed, data$ctd_timed[,5] == contact_types[i])

      ## number of unique places for that contact type - we add one two account
      ## for 'unknown place' (i.e. place = 0 in timeline)
      unq <- length(place_order[[i]])

      ## construct timeline (each row is a case, each column is a day)
      mat <- matrix(-1, nrow = data$N, ncol = t_range)

      ## check p_place
      if(!is.null(data$p_place[[i]])) {

        if(!length(data$p_place[[i]]) %in% c(unq, unq + 1)) {
          stop(sprintf(paste("%d probabilities provided in p_place, but",
                             "%d required (one for each place)"),
                       length(data$p_place[[i]]), unq))
        }
        if(!all.equal(sum(data$p_place[[i]]), 1)) {
          stop("Probabilities in p_place must sum to 1")
        }

      } else {

        ## if not provided assume equal prior probability of all places
        data$p_place[[i]] <- rep(1/unq, unq)

      }

      if(is.null(data$p_trans[[i]])) {

        ## calculate transition matrices for all pairwise combinations of observed places
        comb <- expand.grid(1:unq, 1:unq)
        data$p_trans[[i]] <- matrix(apply(comb, 1, get_trans, data$p_place[[i]]), unq)

      } else if(is.matrix(data$p_trans[[i]])) {

        if(ncol(data$p_trans[[i]]) != unq | nrow(data$p_trans[[i]]) != unq) {
          stop(sprintf(paste("%dx%d matrix of transition probabilities provided, but",
                             "%dx%d required (one row for each place)"),
                       ncol(data$p_trans[[i]]), nrow(data$p_trans[[i]]), unq, unq))
        }

        if(!all(vapply(apply(data$p_trans[[i]], 1, sum), all.equal, TRUE, 1))) {
          stop("Rows of transition matrices must sum to 1")
        }

      } else {

        stop("p_trans must be NULL or a matrix")

      }

      ## update p_place and p_trans to account for 'unobserved' places - these
      ## essentially represent the average of the places we do know of. rmean
      ## is the average probability *from* a given place to any other - it's
      ## therefore just the proportion unobserved divided by the number
      ## unobserved

      ## infer number of unobserved places
      data$N_place_unobserved[i] <- get_n_unobs(unq, data$prop_place_observed[i])

      p_unob = ifelse(data$prop_place_observed[i] < 1,
      (1 - data$prop_place_observed[i])/data$N_place_unobserved[i],
      0)

      rmean <- matrix(c(rep(p_unob, unq), 0), ncol = 1)

      ## cmean is the average probability *to* a given place, from all possible
      ## starting places - it has to be weighted by the prior probability of
      ## starting in any given place - and then has to be weighted by p_observed

      .f <- function(i, p_place, p_trans) {
        p_place <- p_place[-i]
        x <- p_trans[-i,i]
        return(sum(x*p_place))
      }

      cmean <- vapply(seq_along(data$p_place[[i]]), .f,
                      1.0, data$p_place[[i]], data$p_trans[[i]])
      cmean <- cmean*data$prop_place_observed[i]

      ## adjust transition probability to account for unobserved places - the
      ## probability of going to any observed place has to be scaled by
      ## prop_observed
      data$pp_trans[[i]] <- data$p_trans[[i]]*data$prop_place_observed[[i]]

      ## remember that diagonal of transition matrix is 0
      data$pp_trans[[i]] <- rbind(data$pp_trans[[i]],
                                 matrix(cmean, nrow = 1))
      data$pp_trans[[i]] <- cbind(data$pp_trans[[i]], rmean)

      ## also create adjusted transition matrix where the last column represents
      ## the probabilities of transitioning to *any* unobserved place
      tmp <- data$pp_trans[[i]]
      tmp[,ncol(tmp)] <- tmp[,ncol(tmp)]*data$N_place_unobserved[i]
      data$pp_trans_adj[[i]] <- tmp

      ## update prior probability to include unobserved places
      data$pp_place[[i]] <- c(data$p_place[[i]]*data$prop_place_observed[i], p_unob)

      ## also create adjusted prior where the last term is the probability of
      ## being in *any* unobserved place
      tmp <- data$pp_place[[i]]
      tmp[length(tmp)] <- tmp[length(tmp)]*data$N_place_unobserved[i]
      data$pp_place_adj[[i]] <- tmp

      ## fill timeline with places (unknown place = 0)
      for(j in 1:nrow(sub)) {
        ind <- (sub[j,3] + data$C_ind + 1):(sub[j,4] + data$C_ind + 1)
        val <- if(is.na(sub[j, 2])) 0 else match(sub[j, 2], place_order[[i]])
        mat[match(sub[j, 1], data$ids), ind] <- val
      }

      data$ctd_timed_matrix[[i]] <- mat

    }

    ## number of unique places for each contact type
    data$N_place <- as.integer(vapply(data$p_trans, nrow, 1))
    data$N_times <- as.integer(t_range)

    data$has_ctd_timed <- TRUE

  } else {
    data$ctd_timed_matrix <- data$pp_place <- data$pp_place_adj <- list()
    data$pp_trans <- data$pp_trans_adj <- matrix(0, nrow = 0, ncol = 0)
    data$prop_place_observed <-
      data$N_place <-
        data$N_times <-
          data$N_place_unobserved <- numeric()
    data$has_ctd_timed <- FALSE
  }

  ## C++ likelihoods expect this after outbreaker_data (see update_data_with_config)
  if (is.null(data$genetic_model)) {
    data$genetic_model <- "default"
  }

  ## output is a list of checked data
  return(data)

}
