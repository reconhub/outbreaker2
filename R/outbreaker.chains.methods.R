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





#' @export
#' @rdname outbreaker.chains
#' @param what a character string indicating which result to plot
#' @param type a character string indicating the kind of plot to be used: 'series' for time series, 'density' for a kernel density estimation
#' @param burnin the number of iterations to be discarded as burnin
#' @param dens.all a logical indicating if the overal density computed over all runs should be displayed; defaults to TRUE
#' @param col the colors to be used for different runs
#' @param main the title of the plot
#' @param legend a logical indicating if a legend should be added to the plot; defaults to TRUE
#' @param posi a character string indicating the position of the legend; defaults to "bottomleft"
plot.outbreaker.chains <- function(x, what="post", type=c("series","density"), burnin=0, dens.all=TRUE,
                       col=funky(x$n.runs), main=what,
                       legend=TRUE, posi="bottomleft", ...){
    ## HANDLE ARGUMENTS ##
    type <- match.arg(type)
    ## n.runs <- x$n.runs
    n.runs <- 1
    col.ori <- col
    if(!what %in% names(x)) stop(paste(what,"is not a column of x"))
    if(!is.null(col)) col <- rep(col, length = n.runs)
    if(!is.null(lty)) lty <- rep(lty, length = n.runs)
    if(!is.null(lwd)) lwd <- rep(lwd, length = n.runs)
    if(is.null(burnin)){
        burnin <- max(x$burnin, x$find.import.at, x$tune.end)
    }

    ## GET DATA TO PLOT ##
    dat <- cbind(x$step[x$run==1],data.frame(split(x[,what], x$run)))
    names(dat) <- c("step", paste(what, 1:n.runs,sep=""))
    if(!any(dat$step>burnin)) stop("burnin is greater than the number of steps in x")
    dat <- dat[dat$step>burnin,,drop=FALSE]

    ## MAKE PLOT ##
    if(type=="series"){
        matplot(dat$step, dat[,-1,drop=FALSE], type="l", col=col, lty=lty, xlab="MCMC iteration", ylab="value", main=main, ...)
    }

    if(type=="density"){
        ## add general density if needed ##
        temp <- lapply(dat[, -1, drop=FALSE], density)
        if(dens.all){
            temp[[n.runs+1]] <- density(unlist(dat[,-1,drop=FALSE]))
            col <- c(col, "black")
            n.runs <- n.runs+1
        }
        range.x <- range(sapply(temp, function(e) e$x))
        range.y <- range(sapply(temp, function(e) e$y))
        plot(1,type="n", xlim=range.x, ylim=range.y, xlab="value", ylab="density", main=main, ...)
        invisible(lapply(1:n.runs, function(i) lines(temp[[i]], col=col[i], lty=lty[i], lwd=lwd[i])))
    }

    ## ADD LEGEND ##
    if(legend){
        legend(posi, fill=col, title="Runs", leg=1:length(col.ori))
    }
    return(invisible())
} # end plot.outbreaker.chains

