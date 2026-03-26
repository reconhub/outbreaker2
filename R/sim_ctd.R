#' Simulate contact data from a transmission tree
#'
#' Simulates contact data from a transmission tree. The model assumes that all
#' transmission pairs may be observed as contacts and that there is no
#' false-positive reporting of contacts. The probability of contact occurring
#' between a non-transmission pair is given by \code{lambda}. The probability of
#' reporting a contact (transmission pair or not) is given by \code{eps}.
#'
#' @param ttree A data frame or matrix with two columns, each row giving the IDs
#'   (numeric or character) of a transmission pair.
#'
#' @param eps Contact reporting coverage: probability of reporting a contact
#'   (whether or not the pair is a transmission pair).
#'
#' @param lambda Non-infectious contact rate: probability of contact between a
#'   non-transmission pair.
#'
#' @return A data frame with columns \code{i} and \code{j} containing simulated
#'   contact pairs.
#'
#' @author Finlay Campbell (\email{finlaycampbell93@@gmail.com}).
#'
#' @keywords internal
#'
#'
#' @examples
#'
#' ## load data
#' x <- fake_outbreak
#' ttree <- data.frame(from = x$ances, to = seq_along(x$ances))
#'
#' ## simulate contact data with coverage of 80% and 10% probability of
#' ## contact between non-transmission pairs
#' ctd <- outbreaker2:::sim_ctd(ttree, eps = 0.8, lambda = 0.1)
#'
#' ## inspect contact data
#' head(ctd)
sim_ctd <- function(ttree, eps, lambda) {
  ttree <- apply(ttree, 2, as.character)

  if (any(c(eps, lambda) < 0) || any(c(eps, lambda) > 1)) {
    stop("eps and lambda must be probabilities")
  }

  if (ncol(ttree) != 2) {
    stop("ttree must have two columns")
  }

  id <- unique(c(ttree[, 1], ttree[, 2]))
  id <- id[!is.na(id)]

  ## Sort ttree by value or alphabetically, This ensures A:B and B:A are both
  ## recognised when querying the contacts dataframe for transmission pairs
  ttree <- ttree |>
    stats::na.omit() |>
    apply(1, sort, decreasing = FALSE) |>
    t() |>
    as.data.frame(stringsAsFactors = FALSE)

  ttree <- ttree[order(ttree[, 1]), ]
  names(ttree) <- c("V1", "V2")
  if (nrow(ttree) == 0) stop("No transmission observed")

  ## Create a dataframe of all potential contacts
  contacts <- as.data.frame(t(utils::combn(sort(id), 2)),
    stringsAsFactors = FALSE
  )

  ## Create a column of logicals indicating whether a pair represents a
  ## transmission pair or not
  ttree$tPair <- TRUE

  ## Mark non-transmission pairs in the contacts dataframe. The merge function
  ## will mark pairs found in contacts but not in ttree as 'NA'. These are
  ## then converted to FALSE
  contacts <- merge(contacts, ttree, by = c("V1", "V2"), all.x = TRUE)
  contacts$tPair[is.na(contacts$tPair)] <- FALSE

  ## Sample a number of rows given by a binomial distribution
  sampler <- function(x, prob) {
    x[sample(seq_len(nrow(x)), stats::rbinom(1, nrow(x), prob)), 1:3]
  }

  ## Sample transmission pairs with probability eps
  ## Sample non-transmission pairs with probability eps*lambda
  ctd <- rbind(
    sampler(contacts[contacts$tPair, ], eps),
    sampler(contacts[!contacts$tPair, ], eps * lambda)
  )

  ctd <- ctd[, c(1, 2)]

  colnames(ctd) <- c("i", "j")
  rownames(ctd) <- NULL

  return(ctd)
}
