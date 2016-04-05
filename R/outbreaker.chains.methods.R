#' Basic methods for processing outbreaker results
#'
#' Several methods are defined for instances of the class \code{outbreaker.chains}, returned by \code{\link{outbreaker}}, including: \code{print}, \code{plot}
#'
#' @rdname outbreaker.chains
#' @aliases outbreaker.chains
#' @author Thibaut Jombart (\email{thibautjombart@@gmail.com})
#'
#' @param x an \code{outbreaker.chains} object as returned by \code{outbreaker}.
#' @param n.row the number of rows to display in head and tail; defaults to 3.
#' @param n.col the number of columns to display; defaults to 8.
#'
#' @export
#' @importFrom utils head tail
#'
print.outbreaker.chains <- function(x, n.row=3, n.col=8, ...){
    cat("\n\n ///// outbreaker results ///\n")
    cat("\nclass: ", class(x))
    cat("\ndimensions", nrow(x), "rows, ", ncol(x), "columns")

    ## process names of variables not shown
    if(ncol(x) > n.col){
        ori.names <- names(x)
        x <- x[,1:min(n.col, ncol(x))]

        not.shown <- setdiff(ori.names, names(x))

        alpha.txt <- paste(not.shown[range(grep("alpha", not.shown))], collapse=" - ")
        t.inf.txt <- paste(not.shown[range(grep("t.inf", not.shown))], collapse=" - ")

        cat("\nancestries not shown:", alpha.txt)
        cat("\ninfection dates not shown:", t.inf.txt)
    }

    ## heads and tails
    cat("\n\n/// head //\n")
    print(head(as.data.frame(x), n.row))
    cat("\n...")
    cat("\n/// tail //\n")
    print(tail(as.data.frame(x), n.row))
} # end print.outbreaker.chains
