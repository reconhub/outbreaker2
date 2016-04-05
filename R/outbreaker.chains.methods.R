#' Basic methods for processing outbreaker results
#'
#' Several methods are defined for instances of the class \code{outbreaker.chains}, returned by \code{\link{outbreaker}}, including: \code{print}, \code{plot}
#'
#' @rdname outbreaker.chains
#'
#' @aliases outbreaker.chains print.outbreaker.chains plot.outbreaker.chains
#'
#' @author Thibaut Jombart (\email{thibautjombart@@gmail.com})
#'
#' @param x an \code{outbreaker.chains} object as returned by \code{outbreaker}.
#' @param n.row the number of rows to display in head and tail; defaults to 3.
#' @param n.col the number of columns to display; defaults to 8.
#' @param ... further arguments to be passed to other methods
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





#' @rdname outbreaker.chains
#' @param y a character string indicating which result to plot
#' @param type a character string indicating the kind of plot to be used: 'trace' for the MCMC trace, 'hist' for histograms, 'density' for a kernel density estimation
#' @param burnin the number of iterations to be discarded as burnin
## #' @param dens.all a logical indicating if the overal density computed over all runs should be displayed; defaults to TRUE
## #' @param col the colors to be used for different runs
#'
#' @export
#'
#' @importFrom ggplot2 ggplot geom_line geom_point geom_histogram geom_density aes_string labs
#'
plot.outbreaker.chains <- function(x, y="post", type=c("trace", "hist", "density"), burnin=0, ...){

    ## CHECKS ##
    type <- match.arg(type)
    if(!y %in% names(x)) stop(paste(y,"is not a column of x"))


    ## GET DATA TO PLOT ##
    if(burnin > max(x$step)) stop("burnin exceeds the number of steps in x")
    x <- x[x$step>burnin,,drop=FALSE]

    ## MAKE PLOT ##
    if(type=="trace"){
        out <- ggplot(x) + geom_line(aes_string(x="step", y=y)) +
            labs(x="Iteration", y=y, title=paste("trace:",y))
    }

    if(type=="hist"){
        out <- ggplot(x) + geom_histogram(aes_string(x=y)) +
            geom_point(aes_string(x=y, y=0), shape="|", alpha=0.5, size=3) +
        labs(x=y, title=paste("histogram:",y))
    }

    if(type=="density"){
        out <- ggplot(x) + geom_density(aes_string(x=y)) +
            geom_point(aes_string(x=y, y=0), shape="|", alpha=0.5, size=3) +
            labs(x=y, title=paste("density:",y))
    }


    return(out)
} # end plot.outbreaker.chains

