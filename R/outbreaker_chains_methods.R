#' Basic methods for processing outbreaker results
#'
#' Several methods are defined for instances of the class
#' \code{outbreaker_chains}, returned by \code{\link{outbreaker}}, including:
#' \code{print}, \code{plot}
#'
#' @rdname outbreaker_chains
#'
#' @aliases outbreaker_chains print.outbreaker_chains plot.outbreaker_chains
#'   summary.outbreaker_chains
#'
#' @author Thibaut Jombart (\email{thibautjombart@@gmail.com})
#'
#' @param x an \code{outbreaker_chains} object as returned by \code{outbreaker}.
#' @param n_row the number of rows to display in head and tail; defaults to 3.
#' @param n_col the number of columns to display; defaults to 8.
#' @param ... further arguments to be passed to other methods
#'
#' @export
#' @importFrom utils head tail
#'
print.outbreaker_chains <- function(x, n_row = 3, n_col = 8, ...) {
  cat("\n\n ///// outbreaker results ///\n")
  cat("\nclass: ", class(x))
  cat("\ndimensions", nrow(x), "rows, ", ncol(x), "columns")

  ## process names of variables not shown
  if (ncol(x) > n_col) {
    ori_names <- names(x)
    x <- x[, seq_len(min(n_col, ncol(x)))]

    not_shown <- setdiff(ori_names, names(x))

    alpha_txt <- paste(not_shown[range(grep("alpha", not_shown))], collapse=" - ")
    t_inf_txt <- paste(not_shown[range(grep("t_inf", not_shown))], collapse=" - ")
    kappa_txt <- paste(not_shown[range(grep("kappa", not_shown))], collapse=" - ")

    cat("\nancestries not shown:", alpha_txt)
    cat("\ninfection dates not shown:", t_inf_txt)
    cat("\nintermediate generations not shown:", kappa_txt)
  }

  ## heads and tails
  cat("\n\n/// head //\n")
  print(head(as.data.frame(x), n_row))
  cat("\n...")
  cat("\n/// tail //\n")
  print(tail(as.data.frame(x), n_row))
}





#' @rdname outbreaker_chains
#'
#' @param y a character string indicating which result to plot
#'
#' @param type a character string indicating the kind of plot to be used (see details)
#'
#' @param burnin the number of iterations to be discarded as burnin
#'
#' @param min_support a number between 0 and 1 indicating the minimum support of
#' ancestries to be plotted; only used if 'type' is 'network'
#'
## #' @param dens_all a logical indicating if the overal density computed over
## all runs should be displayed; defaults to TRUE #' @param col the colors to be
## used for different runs
#'
#' @export
#'
#' @details 'trace' for the MCMC trace, 'hist' for histograms, 'density' for a
#' kernel density estimation
#'
#' @importFrom ggplot2 ggplot geom_line geom_point geom_histogram geom_density
#' geom_violin aes aes_string coord_flip labs guides scale_size_area
#'
#' @importFrom grDevices xyTable
plot.outbreaker_chains <- function(x, y = "post",
                                   type = c("trace", "hist", "density",
                                            "alpha", "t_inf", "kappa", "network"),
                                   burnin = 0, min_support = 0.1, ...) {

  ## CHECKS ##
  type <- match.arg(type)
  if (!y %in% names(x)) {
    stop(paste(y,"is not a column of x"))
  }

  ## THIS IS JUST TO APPEASE R CMD check

  ## hopefully cran will avoid spurious warnings along the lines of "no
  ## visible binding for global variable" when using ggplot2::aes(...)
  ##
  frequency <- NULL

  ## GET DATA TO PLOT ##
  if (burnin > max(x$step)) {
    stop("burnin exceeds the number of steps in x")
  }
  x <- x[x$step>burnin,,drop = FALSE]

  ## MAKE PLOT ##
  if (type == "trace") {
    out <- ggplot(x) + geom_line(aes_string(x="step", y = y)) +
      labs(x="Iteration", y = y, title = paste("trace:",y))
  }

  if (type == "hist") {
    out <- ggplot(x) + geom_histogram(aes_string(x = y)) +
      geom_point(aes_string(x = y, y = 0), shape="|", alpha = 0.5, size = 3) +
      labs(x = y, title = paste("histogram:",y))
  }

  if (type == "density") {
    out <- ggplot(x) + geom_density(aes_string(x = y)) +
      geom_point(aes_string(x = y, y = 0), shape="|", alpha = 0.5, size = 3) +
      labs(x = y, title = paste("density:",y))
  }

  if (type=="alpha") {
    alpha <- as.matrix(x[,grep("alpha", names(x))])
    colnames(alpha) <- seq_len(ncol(alpha))
    from <- as.vector(alpha)
    to <- as.vector(col(alpha))
    from[is.na(from)] <- 0
    out_dat <- data.frame(xyTable(from,to))
    names(out_dat) <- c("from", "to", "frequency")
    ## Calculate proportion among ancestries
    get.prop <- function(i) {
        ind <- which(out_dat$to == out_dat$to[i])
        out_dat[[3]][i]/sum(out_dat[[3]][ind])
    }
    out_dat[3] <- vapply(seq_along(out_dat[[3]]), get.prop, 1)
    out <- ggplot(out_dat) +
      geom_point(aes(x = factor(to), y = factor(from), size = frequency, color = from)) +
      scale_size_area() +
      guides(colour = FALSE)
  }

  if (type=="t_inf") {
    t_inf <- as.matrix(x[,grep("t_inf", names(x))])
    dates <- as.vector(t_inf)
    cases <- as.vector(col(t_inf))
    out_dat <- data.frame(cases = factor(cases), dates = dates)
    out <- ggplot(out_dat) +
      geom_violin(aes(x = cases, y = dates, fill = cases)) +
      coord_flip() + guides(fill = FALSE) +
      labs(title="infection times")
  }

  if (type=="kappa") {
    kappa <- as.matrix(x[,grep("kappa", names(x))])
    generations <- as.vector(kappa)
    cases <- as.vector(col(kappa))
    to_keep <- !is.na(generations)
    generations <- generations[to_keep]
    cases <- cases[to_keep]
    out_dat <- data.frame(xyTable(generations, cases))
    get.prop <- function(i) {
        ind <- which(out_dat$y == out_dat$y[i])
        out_dat[[3]][i]/sum(out_dat[[3]][ind])
    }
    out_dat[3] <- vapply(seq_along(out_dat[[3]]), get.prop, 1)
    names(out_dat) <- c("generations", "cases", "frequency")
    out <- ggplot(out_dat) +
      geom_point(aes(x = generations, y = cases, size = frequency, color = factor(cases))) +
      scale_size_area() +
      guides(colour = FALSE) +
      labs(title="number of generations between cases", x="number of generations to ancestor")
  }

  if (type=="network") {
    ## extract edge info: ancestries
    alpha <- x[, grep("alpha",names(x)), drop = FALSE]
    from <- unlist(alpha)
    to <- as.vector(col(alpha))
    N <- ncol(alpha)
    edges <- stats::na.omit(data.frame(xyTable(from, to)))
    edges[3] <- edges$number/nrow(alpha)
    names(edges) <- c("from", "to", "value")
    edges <- edges[edges$value > min_support,,drop = FALSE]
    edges$arrows <- "to"
    case_cols <- cases_pal(N)
    edges$color <- case_cols[edges$from]


    ## ## extract edge info: timing
    ## t_inf <- x[, grep("t_inf",names(x)), drop = FALSE]
    ## mean_time <- apply(t_inf, 2, mean)
    ## mean_delay <- mean_time[edges$to] - mean_time[edges$from]
    ## mean_delay[mean_delay<1] <- 1
    ## edges$label <- paste(round(mean_delay), "days")

    ## node info
    find_nodes_size <- function(i) {
      sum(from==i, na.rm = TRUE) / nrow(alpha)
    }
    nodes <- data.frame(id = seq_len(ncol(alpha)),
                        label = seq_len(ncol(alpha)))
    nodes$value <- vapply(nodes$id,
                          find_nodes_size,
                          numeric(1))
    nodes$color <- case_cols
    nodes$shape <- rep("dot", N)

    smry <- summary(x, burnin = burnin)
    is_imported <- is.na(smry$tree$from)
    nodes$shaped[is_imported] <- "star"

    ## generate graph
    out <- visNetwork::visNetwork(nodes = nodes, edges = edges, ...)
    out <- visNetwork::visNodes(out, shadow = list(enabled = TRUE, size = 10),
                                color = list(highlight = "red"))
    out <- visNetwork::visEdges(out, arrows = list(
      to = list(enabled = TRUE, scaleFactor = 0.2)),
      color = list(highlight = "red"))

  }


  return(out)
}




#' @rdname outbreaker_chains
#' @param object an \code{outbreaker_chains} object as returned by \code{outbreaker}.
#' @export
#' @importFrom stats median
summary.outbreaker_chains <- function(object, burnin = 0, ...) {
  ## check burnin ##
  x <- object
  if (burnin > max(x$step)) {
    stop("burnin exceeds the number of steps in object")
  }
  x <- x[x$step>burnin,,drop = FALSE]


  ## make output ##
  out <- list()


  ## summary for $step ##
  interv <- ifelse(nrow(x)>2, diff(tail(x$step, 2)), NA)
  out$step <- c(first = min(x$step),
                last = max(x$step),
                interval = interv,
                n_steps = length(x$step)
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
    as.integer(names(sort(table(x, exclude = NULL), decreasing = TRUE)[1]))
  }
  out$tree$from <- apply(alpha, 2, f1)
  out$tree$to <- seq_len(ncol(alpha))

  ## summary of t_inf ##
  t_inf <- as.matrix(x[,grep("t_inf", names(x))])
  out$tree$time <- apply(t_inf, 2, median)

  ## function to get frequency of most frequent item
  f2 <- function(x) {
    (sort(table(x), decreasing = TRUE)/length(x))[1]
  }
  out$tree$support <- apply(alpha, 2, f2)

  ## summary of kappa ##
  kappa <- as.matrix(x[,grep("kappa", names(x))])
  out$tree$generations <- apply(kappa, 2, median)

  ## shape tree as a data.frame
  out$tree <- as.data.frame(out$tree)
  rownames(out$tree) <- NULL

  return(out)
}
