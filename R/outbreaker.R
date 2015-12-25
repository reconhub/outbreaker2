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
#'
#' @param config a set of settings as returned by \code{\link{outbreaker.config}}
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
outbreaker <- function(data=outbreaker.data(),
                       config=outbreaker.config()){

    ## CHECKS / PROCESS DATA ##
    data <- outbreaker.data(data=data)

    ## CHECK / PROCESS CONFIG ##
    config <- outbreaker.config(data=data, config=config)

    ## INITIALIZE MCMC CHAIN ##
    param <- outbreaker.mcmc.init(data=data, config=config)

    ## CREATE FAST RANDOM VARIABLE GENERATORS ##
    rand <- outbreaker.rand.vec(config)


    ## MCMC ##
    for(i in seq.int(2, config$n.iter, 1)){
        ## move infection dates ##
        if(config$move.t.inf){
            param$current.t.inf <- move.t.inf(data=data, param=param, rand=rand)
        }

        ## move ancestries ##
        if(config$move.ances){
            param$current.ances <- move.ances(data=data, param=param, config=config, rand=rand)
        }

        ## swap ancestries ##
        if(config$move.ances && config$move.t.inf){
            param <- move.swap.ances(data=data, param=param, config=config, rand=rand)
        }

        ## move mu ##
        if(config$move.mu){
            param$current.mu <- move.mu(data=data, param=param, rand=rand)
        }

        ## store outputs if needed
        if((i %% config$sample.every) == 0){
            param <- outbreaker.mcmc.store(chain=param, data=data)
        }

    } # end of the chain


    ## SHAPE RESULTS ##
    param <- outbreaker.mcmc.shape(chain=param, data=data)

    ## RETURN ##
    return(param)
} # end outbreaker

