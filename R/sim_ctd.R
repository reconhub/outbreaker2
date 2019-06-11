#' Simulate contact data from a transmission tree
#'
#' This function simulates contact data from a transmission tree. The
#' model assumes that all transmission pairs have experienced contact, and that
#' there is no false-positive reporting of contacts. The probability of contact
#' occuring between a non-transmision pair is given by the parameter lambda. The
#' probability of reporting a contact (transmission pair or not) is given by the
#' parameters eps.
#'
#' @importFrom magrittr %>%
#'
#' @param tTree a dataframe or matrix of two columns, with each row providing
#'     the ids (numerical or as characters) of a transmission pair
#'
#' @param eps the contact reporting coverage, defined as the probability of
#'     reporting a contact (transmission pair or not)
#'
#' @param lambda the non-infectious contact rate, defined as the probability
#'     of contact between a non-transmission pair.
#'
#' @author Finlay Campbell (\email{f.campbell15@@imperial.ac.uk})
#'
#' @examples
#'
#' ## load data
#' x <- fake_outbreak
#' tTree <- data.frame(from = x$ances, to = seq_along(x$ances))
#'
#' ## simulate contact data with coverage of 80% and 10% probability of
#' ## contact between non-transmission pairs
#' ctd <- outbreaker2:::sim_ctd(tTree, eps = 0.8, lambda = 0.1)
#'
#' ## inspect contact data
#' head(ctd)

sim_ctd <- function(tTree, eps, lambda) {

  tTree <- apply(tTree, 2, as.character)
  
  if(any(c(eps, lambda) < 0) | any(c(eps, lambda) > 1)) {
    stop('eps and lambda must be probabilities')
  }

  if(ncol(tTree) != 2) {
    stop("tTree must have two columns")
  }

  id <- unique(c(tTree[,1], tTree[,2]))
  id <- id[!is.na(id)]
  
  ## Sort tTree by value or alphabetically, This ensures A:B and B:A are both
  ## recognised when querying the contacts dataframe for transmission pairs
  tTree <- tTree %>%
    stats::na.omit() %>%
    apply(1, sort, decreasing = FALSE) %>%
    t() %>%
    as.data.frame(stringsAsFactors = FALSE)

  tTree <- tTree[order(tTree[1]),]
  names(tTree) <- c('V1', 'V2')
  if(nrow(tTree) == 0) stop("No transmission observed")

  ## Create a dataframe of all potential contacts
  contacts <- as.data.frame(t(utils::combn(sort(id), 2)),
                            stringsAsFactors = FALSE)

  ## Create a column of logicals indicating whether a pair represents a
  ## transmission pair or not
  tTree$tPair <- TRUE

  ## Mark non-transmission pairs in the contacts dataframe. The merge function
  ## will mark pairs found in contacts but not in tTree as 'NA'. These are
  ## then converted to FALSE
  contacts <- merge(contacts, tTree, by = c('V1', 'V2'), all.x = TRUE)
  contacts$tPair[is.na(contacts$tPair)] <- FALSE

  ## Sample a number of rows given by a binomial distribution
  sampler <- function(x, prob) {
    x[sample(1:nrow(x), stats::rbinom(1, nrow(x), prob)), 1:3]
  }

  ## Sample transmission pairs with probability eps
  ## Sample non-transmission pairs with probability eps*lambda
  ctd <- rbind(sampler(contacts[contacts$tPair,], eps),
               sampler(contacts[!contacts$tPair,], eps*lambda))

  ctd <- ctd[, c(1, 2)] 
  
  colnames(ctd) <- c('i', 'j')
  rownames(ctd) <- NULL

  return(ctd)
  
}
