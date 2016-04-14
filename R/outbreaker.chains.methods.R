#' Basic methods for processing outbreaker results
#'
#' Several methods are defined for instances of the class \code{outbreaker.chains}, returned by \code{\link{outbreaker}}, including: \code{print}, \code{plot}
#'
#' @rdname outbreaker.chains
#'
#' @aliases outbreaker.chains print.outbreaker.chains plot.outbreaker.chains summary.outbreaker.chains
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
        kappa.txt <- paste(not.shown[range(grep("kappa", not.shown))], collapse=" - ")

        cat("\nancestries not shown:", alpha.txt)
        cat("\ninfection dates not shown:", t.inf.txt)
        cat("\nintermediate generatinons not shown:", kappa.txt)
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
#' @param type a character string indicating the kind of plot to be used (see details)
#' @param burnin the number of iterations to be discarded as burnin
#' @param min.support a number between 0 and 1 indicating the minimum support of ancestries to be plotted; only used if 'type' is 'network'
## #' @param dens.all a logical indicating if the overal density computed over all runs should be displayed; defaults to TRUE
## #' @param col the colors to be used for different runs
#'
#' @export
#'
#' @details
#'  'trace' for the MCMC trace, 'hist' for histograms, 'density' for a kernel density estimation
#'
#' @importFrom ggplot2 ggplot geom_line geom_point geom_histogram geom_density geom_violin aes aes_string coord_flip labs guides
#' @importFrom grDevices xyTable
plot.outbreaker.chains <- function(x, y="post",
                                   type=c("trace", "hist", "density", "alpha", "t.inf", "kappa", "network"),
                                   burnin=0, min.support=0.5, ...){

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

    if(type=="alpha"){
        alpha <- as.matrix(x[,grep("alpha", names(x))])
        colnames(alpha) <- 1:ncol(alpha)
        from <- as.vector(alpha)
        to <- as.vector(col(alpha))
        out.dat <- data.frame(xyTable(from,to))
        out.dat[3] <- out.dat[3]/nrow(x)
        names(out.dat) <- c("from", "to", "frequency")

        out <- ggplot(out.dat) +
            geom_point(aes(x=to, y=from, size=frequency, color=factor(from))) +
                guides(colour=FALSE) +
                    labs(title="ancestries")
    }

    if(type=="t.inf"){
        t.inf <- as.matrix(x[,grep("t.inf", names(x))])
        dates <- as.vector(t.inf)
        cases <- as.vector(col(t.inf))
        out.dat <- data.frame(cases=factor(cases), dates=dates)
        out <- ggplot(out.dat) +
            geom_violin(aes(x=cases, y=dates, fill=cases)) +
                coord_flip() + guides(fill=FALSE) +
                    labs(title="infection times")
    }

    if(type=="kappa"){
        kappa <- as.matrix(x[,grep("kappa", names(x))])
        generations <- as.vector(kappa)
        cases <- as.vector(col(kappa))
        out.dat <- data.frame(xyTable(generations,cases))
        out.dat[3] <- out.dat[3]/nrow(x)
        names(out.dat) <- c("generations", "cases", "frequency")

        out <- ggplot(out.dat) +
            geom_point(aes(x=generations, y=cases, size=frequency, color=factor(cases))) +
                guides(colour=FALSE) +
                    labs(title="number of generations between cases", x="number of generations to ancestor")
    }

    if(type=="network"){
        ## extract edge info: ancestries
        alpha <- x[, grep("alpha",names(x)), drop=FALSE]
        from <- unlist(alpha)
        to <- as.vector(col(alpha))
        edges <- na.omit(data.frame(xyTable(from, to)))
        edges[3] <- edges$number/nrow(alpha)
        names(edges) <- c("from", "to", "value")
        edges <- edges[edges$value > min.support,,drop=FALSE]
        edges$arrows <- "to"

        ## ## extract edge info: timing
        ## t.inf <- x[, grep("t.inf",names(x)), drop=FALSE]
        ## mean.time <- apply(t.inf, 2, mean)
        ## mean.delay <- mean.time[edges$to] - mean.time[edges$from]
        ## mean.delay[mean.delay<1] <- 1
        ## edges$label <- paste(round(mean.delay), "days")

        ## node info
        nodes <- list(id=1:ncol(alpha), label=1:ncol(alpha))
        nodes$value <- sapply(nodes$id, function(i) sum(from==i, na.rm=TRUE))/nrow(alpha)

        ## generate graph
        out <- visNetwork::visNetwork(nodes=nodes, edges=edges, ...)
        out <- visNetwork::visNodes(out, shadow = list(enabled = TRUE, size = 10),
                        color = list(highlight = "red"))
        out <- visNetwork::visEdges(out, arrows = list(
                             to = list(enabled = TRUE, scaleFactor = 0.2)),
                        color = list(highlight = "red"))

    }


    return(out)
} # end plot.outbreaker.chains




#' @rdname outbreaker.chains
#' @param object an \code{outbreaker.chains} object as returned by \code{outbreaker}.
#' @export
#' @importFrom stats median
summary.outbreaker.chains <- function(object, burnin=0, ...){
    ## check burnin ##
    x <- object
    if(burnin > max(x$step)) stop("burnin exceeds the number of steps in object")
    x <- x[x$step>burnin,,drop=FALSE]


    ## make output ##
    out <- list()


    ## summary for $step ##
    interv <- ifelse(nrow(x)>2, diff(tail(x$step, 2)), NA)
    out$step <- c(first = min(x$step),
                  last = max(x$step),
                  interval = interv,
                  n.steps = length(x$step)
                  )


    ## summary of post, like, prior ##
    out$post <- summary(x$post)
    out$like <- summary(x$like)
    out$prior <- summary(x$prior)


    ## summary for mu ##
    out$mu <- summary(x$mu)

    ## summary for pi ##
    out$pi <- summary(x$pi)


    ## summary tree ##
    out$tree <- list()

    ## summary of alpha ##
    alpha <- as.matrix(x[,grep("alpha", names(x))])

    ## function to get most frequent item
    f1 <- function(x) {
        as.integer(names(sort(table(x, exclude=NULL), decreasing=TRUE)[1]))
    }
    out$tree$from <- apply(alpha, 2, f1)
    out$tree$to <- 1:ncol(alpha)

    ## summary of t.inf ##
    t.inf <- as.matrix(x[,grep("t.inf", names(x))])
    out$tree$time <- apply(t.inf, 2, median)

    ## function to get frequency of most frequent item
    f2 <- function(x) {
        (sort(table(x), decreasing=TRUE)/length(x))[1]
    }
    out$tree$support <- apply(alpha, 2, f2)

    ## summary of kappa ##
    kappa <- as.matrix(x[,grep("kappa", names(x))])
    out$tree$generations <- apply(kappa, 2, median)

    ## shape tree as a data.frame
    out$tree <- as.data.frame(out$tree)
    rownames(out$tree) <- NULL

    ## RETURN ##
    return(out)
} # end summary.outbreaker.chains
