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
#' @seealso \code{outbreaker.data} to process input data, and \code{outbreaker.config} to process/set up parameters
#'
#' @param dna the DNA sequences in \code{DNAbin} format (see
#' \code{\link[ape]{read.dna}} in the ape package); this can be imported from a
#' fasta file (extension .fa, .fas, or .fasta) using \code{adegenet}'s function
#' \link[adegenet]{fasta2DNAbin}.
#'
#' @param dates a vector indicating the collection dates, provided either as
#' integer numbers or in a usual date format such as \code{Date} or
#' \code{POSIXct} format. By convention, zero will indicate the oldest date.
#'
#' @param w.dens a vector of numeric values indicating the generation time
#' distribution, reflecting the infectious potential of a case t=0, 1, 2, ...
#' time steps after infection. By convention, w.dens[1]=0, meaning that an
#' newly infected patient cannot be instantaneously infectious. If not
#' standardized, this distribution is rescaled to sum to 1.
#'
#' @param f.dens similar to \code{w.dens}, except that this is the distribution
#' of the colonization time, i.e. time interval during which the pathogen can
#' be sampled from the patient.
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
#' out <- outbreaker(dna=dat$dna, dates=dat$onset, w.dens=w, config=list(n.iter=100, sample.every=10))
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

