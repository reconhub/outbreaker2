\name{Mutation rate estimation}
\alias{get.mu}
\title{Derive mutation rate estimation from outbreak's outputs}
\description{
  The function \code{get.mu} is used to obtain a distribution of the
  mutation rate from outbreaker's ouptput (functions \code{outbreaker}
  and \code{outbreaker.parallel}). The mutation rates used in
  outbreaker's model are expressed per generation of infection,
  which can be problematic to interprete biologically. \code{get.mu}
  derives classical estimates of the mutation rate per unit of time,
  with one value being estimated for each chain of the MCMC. By default,
  the mutation rate is expressed in number of nucleotide changes per
  unit time and per genome. If \code{genome.size} is provided, the
  mutation rate is expressed in number of nucleotide changes per
  unit time and per site.
}
\usage{
get.mu(x, burnin=2e4, genome.size=NULL)
}
\arguments{
  \item{x}{the output of \code{outbreaker} or \code{outbreaker.parallel}.}
  \item{burnin}{an integer indicating the number of steps of the MCMC to
    be discarded as burnin period. Defaults to 20,000.}
  \item{genome.size}{the size of the genome; if not provided, mutation
    rate will be expressed in number of mutations per unit of time and
    per genome.}
}
\value{
  A vector of mutation rates derived from the MCMC.
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk}}
\examples{
## load data
data(fakeOutbreak)
attach(fakeOutbreak)

mu <- get.mu(res, genome.size=ncol(dat$dna))
hist(mu, col="grey",
     main="Inferred distribution of mu",
     xlab="mutations/site/day")
abline(v=1e-4,lty=2, lwd=4, col="royalblue")
mtext(side=3, "Dashed line = actual value")

detach(fakeOutbreak)
}
