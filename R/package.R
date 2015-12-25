#'
#' outbreaker2: disease outbreak reconstruction using epidemiological and genetic data
#'
#' This package re-implements the model introduced by Jombart et al. (PLoS Comput. Biol, 2014) for disease outbreak reconstruction using epidemiological and genetic data. This new version is an alternative to the original package 'outbreaker', relying on a more parsimonious use of C code and implementing new, more flexible sampling approaches.
#'
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @name outbreaker.package
#'
#' @import parallel
#'
#' @importFrom utils packageDescription
#'
#' @importFrom ape dist.dna as.DNAbin
#'
#' @importFrom("stats", "na.omit")
#'
NULL
