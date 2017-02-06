#' Simulate contact tracing data from a simOutbreak object
#'
#' This function takes a simOutbreak object and returns a dataframe of contact pairs. The probabilites of reporting infectious and non-infectious contacts are provided as parameters.
#'
#' @importFrom magrittr %>%
#'
#' @param outbreak a \code{simOutbreak} object
#'
#' @param eps the contact reporting coverage, defined as the probability of
#'     reporting a contact between a transmission pair.
#'
#' @param lambda the non-transmission contact rate, defined as the probability
#'     of contact between a non-transmisison pair.
#'
#' @author Finlay Campbell (\email{f.campbell15@@imperial.ac.uk})
#'
#' @export

sim_CTD <- function(outbreak, eps, lambda) {

    if(outbreak$n < 2) stop("No transmission observed")

    accept_reject <- function(pair, C) {

        is_contact <- C[pair[1], pair[2]]

        if(is_contact) {
            return(stats::runif(1, 0, 1) < eps)
        } else {
            return(stats::runif(1, 0, 1) < lambda * eps)
        }
    }

    ## Extract the transmission tree
    tTree <- cbind(outbreak$ances, outbreak$id)

    ## Define a matrix of all potential contacts
    C <- matrix(FALSE, outbreak$n, outbreak$n)

    ## Mark transmission pairs
    for(i in seq_len(nrow(tTree))) {
        pair <- tTree[i,]
        C[pair[[1]], pair[[2]]] <- C[pair[[2]], pair[[1]]] <- TRUE
    }

    ## Describe all potential contacts
    potent_CTD <- t(utils::combn(outbreak$id, 2))

    accept <- apply(potent_CTD, 1, accept_reject, C)

    CTD <- potent_CTD[accept,]
    colnames(CTD) <- c("i", "j")

    return(CTD)
}
