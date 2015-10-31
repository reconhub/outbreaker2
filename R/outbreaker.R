
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
outbreaker <- function(dates, dna=NULL,
                       w.dens, f.dens=w.dens,
                       config=outbreaker.config()){

    ## CHECKS / PROCESS DATA ##
    data <- outbreaker.data(dates=dates, dna=dna, w.dens=w.dens, f.dens=f.dens)

    ## CHECK / PROCESS CONFIG ##
    config <- outbreaker.config(data=data, config=config)

    ## INITIALIZE MCMC CHAIN ##
    chain <- outbreaker.mcmc.init(data=data, config=config)


    ## MCMC ##
    for(i in 2:config$n.iter){
        ## move infection dates ##
        if(config$move.t.inf){
            chain$current.t.inf <- move.t.inf(data=data, config=config, chain=chain)
        }

        ## move ancestries ##
        if(config$move.ances){
            chain$current.ances <- move.ances(data=data, config=config, chain=chain)
        }

        ## move mu ##
        if(config$move.mu){
            chain$current.mu <- move.mu(data=data, config=config, chain=chain)
        }

        ## store outputs if needed
        if((i %% config$sample.every) == 0){
            chain <- outbreaker.mcmc.store(chain=chain, data=data, config=config)
        }

    } # end of the chain


    ## SHAPE RESULTS ##
    chain <- outbreaker.mcmc.shape(chain=chain, data=data)

    ## RETURN ##
    return(chain)
} # end outbreaker

