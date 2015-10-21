\encoding{UTF-8}
\name{outbreaker}
\alias{outbreaker}
\alias{outbreaker.parallel}
\title{Outbreaker: disease outbreak reconstruction using genetic data}
\description{
  \code{outbreaker} is a tool for the reconstruction of disease
  outbreaks using pathogens genome sequences. It relies on a
  probabilistic model of disease transmission which takes the genetic
  diversity, collection dates, duration of pathogen colonization and
  time interval between cases into account. It is embedded in a Bayesian
  framework which allows to estimate the distributions of parameters of
  interest. It currently allows to estimate:
  \itemize{
    \item transmission trees
    \item dates of infection
    \item missing cases in a chain of transmission
    \item mutation rates
    \item imported cases
    \item (indirectly) effective reproduction numbers
  }

  The function \code{outbreaker} is the basic implementation of the
  model. \code{outbreaker.parallel} allows to run several independent
  MCMC in parallel across different cores / processors of the same
  computer. This requires the base package \code{parallel}.

  The spatial module implemented in outbreaker is currently under
  development. Please contact the author before using it.

  For more resources including tutorials, forums, etc., see:
  \url{http://sites.google.com/site/therepiproject/r-pac/outbreaker}

}
\usage{
outbreaker(dna=NULL, dates, idx.dna=NULL, mut.model=1, spa.model=0,
           w.dens, w.trunc=length(w.dens), f.dens=w.dens,
           f.trunc=length(f.dens), dist.mat=NULL,
           init.tree=c("seqTrack","random","star"), init.kappa=NULL,
           init.mu1=NULL, init.mu2=init.mu1, init.spa1=NULL, n.iter=1e5,
           sample.every=500, tune.every=500, burnin=2e4,
           import.method=c("genetic","full","none"), find.import.n=50,
           pi.prior1=10, pi.prior2=1, spa1.prior=1, move.mut=TRUE,
           move.ances=TRUE, move.kappa=TRUE, move.Tinf=TRUE,
           move.pi=TRUE, move.spa=TRUE, outlier.threshold
           = 5, max.kappa=10, quiet=TRUE, res.file.name="chains.txt",
           tune.file.name="tuning.txt", seed=NULL)


outbreaker.parallel(n.runs, parallel=TRUE, n.cores=NULL,
           dna=NULL, dates, idx.dna=NULL, mut.model=1, spa.model=0,
           w.dens, w.trunc=length(w.dens), f.dens=w.dens,
           f.trunc=length(f.dens), dist.mat=NULL,
           init.tree=c("seqTrack","random","star"), init.kappa=NULL,
           init.mu1=NULL, init.mu2=init.mu1, init.spa1=NULL, n.iter=1e5,
           sample.every=500, tune.every=500, burnin=2e4,
           import.method=c("genetic","full","none"), find.import.n=50,
           pi.prior1=10, pi.prior2=1, spa1.prior=1, move.mut=TRUE,
           move.ances=TRUE, move.kappa=TRUE, move.Tinf=TRUE,
           move.pi=TRUE, move.spa=TRUE, outlier.threshold
           = 5, max.kappa=10, quiet=TRUE, res.file.name="chains.txt",
           tune.file.name="tuning.txt", seed=NULL)
}
\arguments{
  \item{dna}{the DNA sequences in \code{DNAbin} format (see
    \code{\link[ape]{read.dna}} in the ape package); this can be
    imported from a fasta file (extension .fa, .fas, or .fasta) using
    \code{adegenet}'s function \link[adegenet]{fasta2DNAbin}.}

  \item{dates}{a vector indicating the collection dates, provided either
    as integer numbers or in a usual date format such as \code{Date} or
    \code{POSIXct} format. By convention, zero will indicate the oldest
    date.}
  
  \item{idx.dna}{an optional integer vector indicating to which case
    each dna sequence in \code{dna} corresponds. Not required if each
    case has a sequence, and the order of the sequences matches that of
    the cases.}

  \item{mut.model}{an integer indicating the mutational model to be
    used; 1: one single mutation rate; 2: two rates,
    transitions (mu1) / transversions (mu2).}
  
  \item{spa.model}{an integer indicating the spatial model to be
    used. 0: no spatial model (default). 1: exponential kernel (under
    development).}
 
  \item{w.dens}{a vector of numeric values indicating the generation
    time distribution, reflecting the infectious potential of a case
    t=0, 1, 2, ... time steps after infection. By convention,
    w.dens[1]=0, meaning that an newly infected patient cannot be
    instantaneously infectious. If not standardized, this distribution
    is rescaled to sum to 1.}
  
  \item{w.trunc}{an integer indicating after which time step the
    distribution \code{w.dens} should be truncated to zero (effectively,
    the maximum duration of the infectious period).}
  
  \item{f.dens}{similar to \code{w.dens}, except that this is the
    distribution of the colonization time, i.e. time interval during
    which the pathogen can be sampled from the patient.}
  
  \item{f.trunc}{an integer indicating after which time step the
    distribution \code{f.dens} should be truncated to zero (effectively,
    the maximum duration of colonization).}

  \item{dist.mat}{a matrix of pairwise spatial distances between the
    cases.}

  %% \item{locations}{a factor indicating the location of individuals used
  %%   to defined local transmissions in the stratified dispersal model}
  
  \item{init.tree}{the tree used to initialize the MCMC. Can be either a
    character string indicating how this tree should be computed, or a
    vector of integers corresponding to the tree itself, where the i-th
    value corresponds to the index of the ancestor of 'i' (i.e.,
    \code{init.tree[i]} is the ancestor of case \code{i}). Accepted
    character strings are "seqTrack" (uses seqTrack output as initialize
    tree), "random" (ancestor randomly selected from preceding cases),
    and "star" (all cases coalesce to the first case). Note that for
    SeqTrack, all cases should have been sequenced.}
  
  \item{init.kappa}{as \code{init.tree}, but values indicate the number
    of generations between each case and its most recent sampled
    ancestor.}

  \item{n.iter}{an integer indicating the number of iterations in the
    MCMC, including the burnin period; defaults to \code{100,000}.}
  
  \item{sample.every}{an integer indicating the frequency at which to
    sample from the MCMC, defaulting to 500 (i.e., output to file every
    500 iterations).}
  
  \item{tune.every}{an integer indicating the frequency at which
    proposal distributions are tuned, defaulting to 500 (i.e., tune
    proposal distribution every 500 iterations).}
  
  \item{burnin}{an integer indicating the number of iterations for the
    burnin period, after which the chains are supposed to have mixed;
    estimated values of parameter are only relevant after the burnin
    period. Used only when imported cases are automatically detected.}

  \item{import.method}{a character string indicating which method to use
    for detecting imported cases; available choices are 'gen' (based on
    genetic likelihood), 'full' (based on full likelihood), and 'none' (no
    imported case detection).}
  
  \item{find.import.n}{an integer indicating how many chains should be
    used to determine imported cases; note that this corresponds to
    chains that are output after the burnin, so that a total of (burnin
    + output.every*find.import.n) chains will be used in the prior run
    to determine imported cases. Defaults to \code{50}.}
  
  \item{pi.prior1, pi.prior2}{two numeric values being the parameters of
    the Beta distribution used as a prior for \eqn{\pi}. This prior
    is Beta(10,1) by default, indicating that a majority of cases are
    likely to have been observed. Use Beta(1,1) for a flat prior.}

  %% \item{phi.param1, phi.param2}{two numeric values being the parameters of
  %%   the Beta distribution used as a prior for \eqn{\phi}. This prior
  %%   is Beta(5,1) by default, indicating that a majority of cases are
  %%   likely to have been observed. Use Beta(1,1) for a flat prior.}

  \item{init.mu1,init.mu2}{initial values for the mutation rates (mu1:
    transitions; mu2: transversions).}

  \item{init.spa1}{initial values of the spatial parameter.}

  \item{spa1.prior}{parameters of the prior distribution for
    the spatial parameters. In the spatial model 1, \code{spa1.prior} is
    the mean of an exponential distribution.}

  \item{move.mut,move.pi,move.spa}{logicals indicating whether
    the named items should be estimated ('moved' in the MCMC), or not,
    all defaulting to TRUE. \code{move.mut} handles both mutation
    rates.}

  \item{move.ances,move.kappa,move.Tinf}{vectors of logicals of length
    'n' indicating for which cases different components should be
    moved during the MCMC.}
  
  \item{outlier.threshold}{a numeric value indicating the threshold for
    detecting low likelihood values corresponding to imported
    cases. Outliers have a likelihood \code{outlier.threshold} smaller
    than the average.}

  \item{max.kappa}{an integer indicating the maximum number of
    generations between a case and its most recent sampled ancestor;
    defaults to 10.}
  
  \item{quiet}{a logical indicating whether messages should be displayed
    on the screen.}
  
  \item{res.file.name}{a character string indicating the name of the
    file used to store MCMC outputs.}
  
  \item{tune.file.name}{a character string indicating the name of the
    file used to store MCMC tuning outputs.}
  
  \item{seed}{an integer used to set the random seed of the C
    procedures.}  \item{n.runs}{an integer indicating the number of
    independent chains to run, either in parallel (if \code{parallel} is
    used), or serially (otherwise).}
  
  \item{parallel}{a logical indicating whether the package
    \code{parallel} should be used to run parallelized computations; by
    default, it is used if available.}
  
  \item{n.cores}{an integer indicating the number of cores to be used
    for parallelized computations; if NULL (default value), then up to 6
    cores are used, depending on availability.}
}
\value{
  Both procedures return a list with the following components:
  \itemize{
    \item chains: a data.frame containing MCMC outputs (which
    are also stored in the file indicated in \code{res.file.name}).

    \item collec.dates: (data) the collection dates.

    \item w: (data) the generation time distribution (argument
    \code{w.dens})

    \item f: (data) the distribution of the time to collection (argument
    \code{f.dens})

    \item D: a matrix of genetic distances (in number of mutations)
    between all pairs of sequences.

    \item idx.dna: (data) the index of the case each dna sequence
    corresponds to

    \item tune.end: an integer indicating at which iteration the
    proposal auto-tuning procedures all stopped.

    \item find.import: a logical indicating if imported cases were to be
    automatically detected.

    \item burnin: an integer indicating the pre-defined burnin, used
    when detecting imported cases.

    \item find.import.at: an integer indicating at which iteration of
    the preliminary MCMC imported cases were detected.

    \item n.runs: the number of independent runs used.

    \item call: the matched call.  }
}
\references{
 Jombart T, Cori A, Didelot X, Cauchemez S, Fraser C and Ferguson N
(accepted).  Bayesian reconstruction of disease outbreaks by combining
epidemiologic and genomic data. PLoS Computational Biology.

}
\author{
  Thibaut Jombart (\email{t.jombart@imperial.ac.uk})
}
\seealso{
  \itemize{
    \item \link{plotChains} to visualize MCMC chains.

    \item \link{transGraph} and \link{get.tTree} to represent
    transmission trees.
    
    \item \link{get.R} and \link{get.Rt} to get reproduction numbers
    distributions.
    
    \item \link{get.incid} to get estimates of incidence.

    \item \link{get.mu} to get the mutation rate distribution.
    
    \item \link{simOutbreak} to simulate outbreaks.
    
    \item \link{selectChains} to select chains from parallel runs which
    converged towards different posterior modes.

    \item \link{fakeOutbreak}, a toy dataset used to illustrate the
    method.

    \item For more resources including tutorials, forums, etc., see:
    \url{http://sites.google.com/site/therepiproject/r-pac/outbreaker}

  }
}
\examples{

## EXAMPLE USING TOYOUTBREAK ##
## LOAD DATA, SET RANDOM SEED
data(fakeOutbreak)
attach(fakeOutbreak)

## VISUALIZE DYNAMICS
matplot(dat$dynam, type="o", pch=20, lty=1,
   main="Outbreak dynamics", xlim=c(0,28))
legend("topright", legend=c("S","I","R"), lty=1, col=1:3)

## VISUALIZE TRANSMISSION TREE
plot(dat, annot="dist", main="Data - transmission tree")
mtext(side=3, "arrow annotations are numbers of mutations")


\dontrun{
## RUN OUTBREAKER - PARALLEL VERSION
## (takes < 1 min))
set.seed(1)
res <-  outbreaker.parallel(n.runs=4, dna=dat$dna,
   dates=collecDates,w.dens=w, n.iter=5e4)
}


## ASSESS CONVERGENCE OF CHAINS
plotChains(res)
plotChains(res, burnin=2e4)

## REPRESENT POSTERIOR ANCESTRIES
transGraph(res, annot="", main="Posterior ancestries", thres=.01)

## GET CONSENSUS ANCESTRIES
tre <- get.tTree(res)
plot(tre, annot="", main="Consensus ancestries")

## SHOW DISCREPANCIES
col <- rep("lightgrey", 30)
col[which(dat$ances != tre$ances)] <- "pink"
plot(tre, annot="", vertex.color=col, main="Consensus ancestries")
mtext(side=3, text="cases with erroneous ancestries in pink")

## GET EFFECTIVE REPRODUCTION OVER TIME
get.Rt(res)

## GET INDIVIDUAL EFFECTIVE REPRODUCTION
head(get.R(res))
boxplot(get.R(res), col="grey", xlab="Case",
        ylab="Effective reproduction number")

## GET MUTATION RATE PER TIME UNIT
## per genome
head(get.mu(res))

## per nucleotide
mu <- get.mu(res, genome.size=1e4)
head(mu)

summary(mu)
hist(mu, border="lightgrey", col="grey", xlab="Mutation per day and nucleotide",
     main="Posterior distribution of mutation rate")

detach(fakeOutbreak)


}