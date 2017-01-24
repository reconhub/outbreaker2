#' @rdname outbreaker_config
#' @export
#' 
#' @aliases print.outbreaker_config
#'
#' @param x an \code{outbreaker_config} object as returned by \code{create_config}.
#' @param ... further arguments to be passed to other methods
#' 
print.outbreaker_config <- function(x, ...) {
    cat("\n\n ///// outbreaker settings ///\n")
    cat("\nclass:", class(x))
    cat("\nnumber of items:", length(x))

    cat("\n\n/// initialisation //\n")
    to_print <- grep("init", names(x))
    print(noquote(unlist(x[to_print])))
    x <- x[-to_print]

    cat("\n/// movements //\n")
    to_print <- unlist(sapply(c("move", "^sd"), grep, names(x)))
    print(noquote(sapply(x[to_print], as.character)))
    x <- x[-to_print]

    cat("\n/// chains //\n")
    to_print <- c("n_iter", "sample_every")
    print(noquote(sapply(x[to_print], as.character)))
    x <- x[-match(to_print, names(x))]

    cat("\n/// priors //\n")
    to_print <- grep("prior", names(x))
    print(noquote(unlist(x[to_print])))
    x <- x[-to_print]

    cat("\n/// imported cases //\n")
    to_print <- unlist(sapply(c("import", "threshold"), grep, names(x)))
    print(noquote(sapply(x[to_print], as.character)))
    x <- x[-to_print]

    if(length(x)>0){
        cat("\n/// other settings //\n")
        print(noquote(sapply(x, as.character)))
    }

    return(invisible(NULL))
}
