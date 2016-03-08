#' Outbreaker2: disease outbreak reconstruction using epidemiological and genetic data
#'
#'
#' @export
#'
#' @aliases outbreaker
#'
#' @rdname outbreaker
#'
#' @author Thibaut Jombart (\email{t.jombart@@imperial.ac.uk})
#'
#' @references Jombart T, Cori A, Didelot X, Cauchemez S, Fraser C and Ferguson
#' N (2014).  Bayesian reconstruction of disease outbreaks by combining
#' epidemiologic and genomic data. PLoS Computational Biology.
#'
#' @seealso \code{\link{outbreaker.data}} to process input data, and \code{\link{outbreaker.config}} to process/set up parameters
#'
#' @param data a list of named items containing input data as returned by \code{\link{outbreaker.data}}
#' @param config a set of settings as returned by \code{\link{outbreaker.config}}
#' @param loglike a set of log-likelihood functions as returned by \code{\link{outbreaker.create.loglike}}
#' @param priors a set of prior functions as returned by \code{\link{outbreaker.create.priors}}
#' @param posteriors a set of posterior functions as returned by \code{\link{outbreaker.create.posteriors}}
#' @param moves a set of movement functions stored in a named list as returned by \code{\link{outbreaker.create.moves}}
#'
#' @seealso \code{\link{outbreaker.config}} to see default parameters / set parameters and \code{\link{outbreaker.data}} to process input data
#'
#' @importFrom stats dexp rnorm runif
#' @importFrom coda mcmc
#'
#' @examples
#'
#' ## get data
#' data(fakeOutbreak)
#' dat <- fakeOutbreak$dat
#' w <- fakeOutbreak$w
#'
#' ## run outbreaker
#' out <- outbreaker(data=list(dna=dat$dna, dates=dat$onset, w.dens=w),
#' config=list(n.iter=100, sample.every=10))
#'
#'
outbreaker <- function(data = outbreaker.data(),
                       config = outbreaker.config(),
                       loglike = outbreaker.create.loglike(),
                       priors = outbreaker.create.priors(),
                       posteriors = outbreaker.create.posteriors(),
                       moves = outbreaker.create.moves()){

    ## CHECKS / PROCESS DATA ##
    data <- outbreaker.data(data=data)

    ## CHECK / PROCESS CONFIG ##
    config <- outbreaker.config(data=data, config=config)

    ## INITIALIZE MCMC CHAIN ##
    param <- outbreaker.mcmc.init(data=data, config=config)

    ## CREATE FAST RANDOM VARIABLE GENERATORS ##
    rand <- outbreaker.rand.vec(config)

    ## MCMC ##
    ## add likelihood, prior and posterior functions to movement context
    add.to.context(environment(), c(loglike=loglike, priors=priors, posteriors=posteriors))

    ## perform mcmc
    param <- outbreaker.move(data=data, config=config, param=param, rand=rand)

    ## SHAPE RESULTS ##
    param <- outbreaker.mcmc.shape(param=param, data=data)

    ## RETURN ##
    return(param)
} # end outbreaker

