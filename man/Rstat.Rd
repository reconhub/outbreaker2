\name{reproduction numbers}
\alias{get.R}
\alias{get.Rt}
\alias{get.incid}
\title{Derive reproduction numbers from outbreak's outputs}
\description{
  These functions are used to compute reproduction numbers and derive
  incidence curves from outbreaker's ouptput (functions
  \code{outbreaker} and \code{outbreaker.parallel}). They all rely on
  the entire outbreak having been sampled.
  \itemize{
    \item \code{get.R} derive distributions of individual effective reproduction numbers.
    \item \code{get.Rt} derives effective reproduction numbers averaged
    for each time step.
    \item \code{get.incid} derives incidence curves for each time step.
  }
}
\usage{
get.Rt(x, burnin=2e4, plot=TRUE, type=c("boxplot", "lines"), lines=FALSE,
       fill.col="gold", lines.col=transp("grey"), \dots)

get.R(x, burnin=2e4, \dots)

get.incid(x, burnin=2e4, plot=TRUE, type=c("boxplot", "lines"),
          lines=FALSE, fill.col="gold", lines.col=transp("grey"), \dots)

}
\arguments{
  \item{x}{the output of \code{outbreaker} or \code{outbreaker.parallel}.}
  \item{burnin}{an integer indicating the number of steps of the MCMC to
    be discarded as burnin period. Defaults to 20,000.}
  \item{plot}{a logical indicating whether a plot should be displayed.}
  \item{type}{a character indicating the type of plot to be used.}
  \item{lines}{a logical indicating whether individual lines should be
    added to the plot.}
  \item{fill.col}{the color to be used for the boxplot.}
  \item{lines.col}{the color to be used to the lines.}
  \item{\dots}{further arguments to be passed to other functions.}
}
\value{
  These functions return a \code{data.frame} containing the plotted information.
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk}}
\examples{
## load data
data(fakeOutbreak)
attach(fakeOutbreak)

## individual R
barplot(table(get.R(res)), main="Individual effective reproduction numbers")

## R(t)
get.Rt(res)

## incidence
get.incid(res)

detach(fakeOutbreak)

}
