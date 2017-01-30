#' outbreaker2: a platform for disease outbreak reconstruction
#'
#' This package provides a statistical platform for reconstructing transmission
#' trees ('who infects whom') in densely sampled disease outbreaks. It
#' reimplements, as a particular case, the model of 'outbreaker' (see
#' references). 'outbreaker2' extends and replaces 'outbreaker'. \cr
#'
#' The emphasis of 'outbreaker2' is on modularity, which enables customisation
#' of priors, likelihoods and even movements of parameters and augmented data by
#' the user. This the dedicated vignette on this topic
#' \code{vignette("outbreaker2_custom")}.\cr
#'
#' The main functions of the package are:
#' \itemize{
#'
#' \item \code{\link{outbreaker}}: main function to run analyses
#'
#' \item \code{\link{outbreaker_data}}: function to process input data
#'
#' \item \code{\link{create_config}}: function to create default and customise
#' configuration settings
#'
#' \item \code{\link{custom_priors}}: function to specify customised prior
#' functions
#'
#' \item \code{\link{custom_likelihoods}}: function to specify customised likelihoods
#' functions
#'
#' \item \code{\link{custom_moves}}: function to create default and customise movement
#' functions
#' 
#' }
#'
#' 
#' This package is in testing state. Please email the authors if you plan on
#' using it.
#'
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @name outbreaker_package
#'
#' @importFrom utils packageDescription
#'
NULL
