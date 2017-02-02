#' Small simulated outbreak
#'
#' This dataset is a small (30 cases) simulated outbreak originally used to
#' illustrate \code{outbreaker}, and used for the same purposes in
#' \code{outbreaker2}.  This list contains the following:
#'
#' \itemize{
#'
#' \item \code{$onset}: A vector of integers representing dates of onset.
#'
#' \item \code{$sample}: A vector of integers representing the dates of isolation.
#'
#' \item \code{$dna}: A DNAbin matrix of pathogen genome sequences.
#'
#' \item \code{$w}: A numeric vector giving the empirical distribution of the
#' generation time.
#'
#' \item \code{$ances}: A vector of integers indicating, for each case, the true
#' infectors. \code{NA} represents imported cases.
#'
#' }
#'
#' @aliases fake_outbreak
#' @docType data
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#' @keywords datasets
#'
#' @examples
#' names(fake_outbreak)
#' fake_outbreak
#' 
"fake_outbreak"


