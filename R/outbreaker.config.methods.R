#' Basic methods for processing outbreaker config
#'
#' Several methods are defined for instances of the class \code{outbreaker.config}, returned by \code{\link{outbreaker.config}}, including: \code{print}
#'
#' @rdname outbreaker.config.methods
#'
#' @aliases print.outbreaker.config
#'
#' @author Thibaut Jombart (\email{thibautjombart@@gmail.com})
#'
#' @param x an \code{outbreaker.config} object as returned by \code{outbreaker.config}.
#' @param ... further arguments to be passed to other methods
#'
#' @export
#'
print.outbreaker.config <- function(x, ...) {
    cat("\n\n ///// outbreaker settings ///\n")
    cat("\nclass:", class(x))
    cat("\nnumber of items:", length(x))

    cat("\n\n/// initialisation //\n")
    to.print <- grep("init", names(x))
    print(noquote(unlist(x[to.print])))
    x <- x[-to.print]

    cat("\n/// movements //\n")
    to.print <- unlist(sapply(c("move", "^sd"), grep, names(x)))
    print(noquote(sapply(x[to.print], as.character)))
    x <- x[-to.print]

    cat("\n/// chains //\n")
    to.print <- c("n.iter", "sample.every")
    print(noquote(sapply(x[to.print], as.character)))
    x <- x[-match(to.print, names(x))]

    cat("\n/// priors //\n")
    to.print <- grep("prior", names(x))
    print(noquote(unlist(x[to.print])))
    x <- x[-to.print]

    cat("\n/// imported cases //\n")
    to.print <- unlist(sapply(c("import", "threshold"), grep, names(x)))
    print(noquote(sapply(x[to.print], as.character)))
    x <- x[-to.print]

    if(length(x)>0){
        cat("\n/// other settings //\n")
        print(noquote(sapply(x, as.character)))
    }

    return(invisible(NULL))
}
