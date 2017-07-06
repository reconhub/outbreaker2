#' @importFrom Rcpp evalCpp
#' @useDynLib outbreaker2, .registration = TRUE
#'
.onAttach <- function(libname, pkgname) {
    ## pkg_version <- packageDescription("outbreaker2", fields = "Version")
    ## startup_txt <- paste("\n   === outbreaker2", pkg_version, "is loaded ===\n")
    ## startup_txt <- paste(startup_txt, "Information & Documentation -> check the R-epi project: \nhttp://sites_google_com/site/therepiproject\n", sep="\n")
    ## startup_txt <- paste(startup_txt, "Questions -> ask the R-epi forum: \nr-epi@googlegroups_com", sep="\n")
    ## packageStartupMessage(startup_txt)
}
