#' Simulate contact tracing data from a transmission tree.
#'
#' This function simulates contact tracing data from a transmission tree. The
#' probability of a contact occuring between a transmision pair is assumed to be
#' 1. The probability of contact occuring between a non-transmision pair is
#' given by the parameter lambda. The probability of reporting a contact
#' (transmission pair or not) is given by the parameters eps.
#'
#' @importFrom magrittr %>%
#'
#' @param tTree a dataframe or matrix of two columns, with each row providing
#'     the ids (numerical or as characters) of a transmission pairs
#'
#' @param eps the contact reporting coverage, defined as the probability of
#'     reporting a contact between a transmission pair.
#'
#' @param lambda the non-transmission contact rate, defined as the probability
#'     of contact between a non-transmisison pair.
#'
#' @author Finlay Campbell (\email{f.campbell15@@imperial.ac.uk})

sim_ctd <- function(tTree, eps, lambda) {

    if(any(c(eps, lambda) < 0) | any(c(eps, lambda) > 1)) {
        stop('eps and lambda must be probabilities')
    }

    id <- unique(c(tTree[,1], tTree[,2]))
    id <- id[!is.na(id)]
    
    ## Sort tTree by value or alphabetically, This ensures A:B and B:A are both
    ## recognised when querying the contacts dataframe for transmission pairs
    tTree <- tTree %>%
        stats::na.omit() %>%
        apply(1, sort, FALSE) %>%
        t() %>%
        as.data.frame(stringsAsFactors = FALSE)

    tTree <- tTree[order(tTree[1]),]
    names(tTree) <- c('V1', 'V2')
    if(nrow(tTree) == 0) stop("No transmission observed")

    ## Create a dataframe of all potential contacts
    contacts <- as.data.frame(t(utils::combn(id, 2)))

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

    rownames(ctd) <- NULL

    return(ctd)
}
