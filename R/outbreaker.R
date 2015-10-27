
#' Outbreaker: disease outbreak reconstruction using genetic data
#'
#' \code{outbreaker} is a tool for the reconstruction of disease outbreaks
#' using pathogens genome sequences. It relies on a probabilistic model of
#' disease transmission which takes the genetic diversity, collection dates,
#' duration of pathogen colonization and time interval between cases into
#' account. It is embedded in a Bayesian framework which allows to estimate the
#' distributions of parameters of interest. It currently allows to estimate:
#' \itemize{ \item transmission trees \item dates of infection \item missing
#' cases in a chain of transmission \item mutation rates \item imported cases
#' \item (indirectly) effective reproduction numbers }
#'
#' The function \code{outbreaker} is the basic implementation of the model.
#' \code{outbreaker.parallel} allows to run several independent MCMC in
#' parallel across different cores / processors of the same computer. This
#' requires the base package \code{parallel}.
#'
#' The spatial module implemented in outbreaker is currently under development.
#' Please contact the author before using it.
#'
#' For more resources including tutorials, forums, etc., see:
#' \url{http://sites.google.com/site/therepiproject/r-pac/outbreaker}
#'
#' @export
#'
#' @aliases outbreaker outbreaker.parallel
#'
#' @rdname outbreaker
#'
#' @param dna the DNA sequences in \code{DNAbin} format (see
#' \code{\link[ape]{read.dna}} in the ape package); this can be imported from a
#' fasta file (extension .fa, .fas, or .fasta) using \code{adegenet}'s function
#' \link[adegenet]{fasta2DNAbin}.
#' @param dates a vector indicating the collection dates, provided either as
#' integer numbers or in a usual date format such as \code{Date} or
#' \code{POSIXct} format. By convention, zero will indicate the oldest date.
#' @param idx.dna an optional integer vector indicating to which case each dna
#' sequence in \code{dna} corresponds. Not required if each case has a
#' sequence, and the order of the sequences matches that of the cases.
#' @param mut.model an integer indicating the mutational model to be used; 1:
#' one single mutation rate; 2: two rates, transitions (mu1) / transversions
#' (mu2).
#' @param spa.model an integer indicating the spatial model to be used. 0: no
#' spatial model (default). 1: exponential kernel (under development).
#' @param w.dens a vector of numeric values indicating the generation time
#' distribution, reflecting the infectious potential of a case t=0, 1, 2, ...
#' time steps after infection. By convention, w.dens[1]=0, meaning that an
#' newly infected patient cannot be instantaneously infectious. If not
#' standardized, this distribution is rescaled to sum to 1.
#' @param f.dens similar to \code{w.dens}, except that this is the distribution
#' of the colonization time, i.e. time interval during which the pathogen can
#' be sampled from the patient.
#' @param dist.mat a matrix of pairwise spatial distances between the cases.
#' @param init.tree the tree used to initialize the MCMC. Can be either a
#' character string indicating how this tree should be computed, or a vector of
#' integers corresponding to the tree itself, where the i-th value corresponds
#' to the index of the ancestor of 'i' (i.e., \code{init.tree[i]} is the
#' ancestor of case \code{i}). Accepted character strings are "seqTrack" (uses
#' seqTrack output as initialize tree), "random" (ancestor randomly selected
#' from preceding cases), and "star" (all cases coalesce to the first case).
#' Note that for SeqTrack, all cases should have been sequenced.
#' @param init.kappa as \code{init.tree}, but values indicate the number of
#' generations between each case and its most recent sampled ancestor.
#' @param n.iter an integer indicating the number of iterations in the MCMC,
#' including the burnin period; defaults to \code{100,000}.
#' @param sample.every an integer indicating the frequency at which to sample
#' from the MCMC, defaulting to 500 (i.e., output to file every 500
#' iterations).
#' @param tune.every an integer indicating the frequency at which proposal
#' distributions are tuned, defaulting to 500 (i.e., tune proposal distribution
#' every 500 iterations).
#' @param burnin an integer indicating the number of iterations for the burnin
#' period, after which the chains are supposed to have mixed; estimated values
#' of parameter are only relevant after the burnin period. Used only when
#' imported cases are automatically detected.
#' @param import.method a character string indicating which method to use for
#' detecting imported cases; available choices are 'gen' (based on genetic
#' likelihood), 'full' (based on full likelihood), and 'none' (no imported case
#' detection).
#' @param find.import.n an integer indicating how many chains should be used to
#' determine imported cases; note that this corresponds to chains that are
#' output after the burnin, so that a total of (burnin +
#' output.every*find.import.n) chains will be used in the prior run to
#' determine imported cases. Defaults to \code{50}.
#' @param pi.prior1,pi.prior2 two numeric values being the parameters of the
#' Beta distribution used as a prior for \eqn{\pi}. This prior is Beta(10,1) by
#' default, indicating that a majority of cases are likely to have been
#' observed. Use Beta(1,1) for a flat prior.
#' @param init.mu1,init.mu2 initial values for the mutation rates (mu1:
#' transitions; mu2: transversions).
#' @param init.spa1 initial values of the spatial parameter.
#' @param spa1.prior parameters of the prior distribution for the spatial
#' parameters. In the spatial model 1, \code{spa1.prior} is the mean of an
#' exponential distribution.
#' @param move.mut,move.pi,move.spa logicals indicating whether the named items
#' should be estimated ('moved' in the MCMC), or not, all defaulting to TRUE.
#' \code{move.mut} handles both mutation rates.
#' @param move.ances,move.kappa,move.Tinf vectors of logicals of length 'n'
#' indicating for which cases different components should be moved during the
#' MCMC.
#' @param outlier.threshold a numeric value indicating the threshold for
#' detecting low likelihood values corresponding to imported cases. Outliers
#' have a likelihood \code{outlier.threshold} smaller than the average.
#' @param max.kappa an integer indicating the maximum number of generations
#' between a case and its most recent sampled ancestor; defaults to 10.
#' @param quiet a logical indicating whether messages should be displayed on
#' the screen.
#' @param res.file.name a character string indicating the name of the file used
#' to store MCMC outputs.
#' @param tune.file.name a character string indicating the name of the file
#' used to store MCMC tuning outputs.
#' @param seed an integer used to set the random seed of the C procedures.
#' @param n.runs an integer indicating the number of independent chains to run,
#' either in parallel (if \code{parallel} is used), or serially (otherwise).
#' @param parallel a logical indicating whether the package \code{parallel}
#' should be used to run parallelized computations; by default, it is used if
#' available.
#' @param n.cores an integer indicating the number of cores to be used for
#' parallelized computations; if NULL (default value), then up to 6 cores are
#' used, depending on availability.
#' @return Both procedures return a list with the following components:
#' \itemize{ \item chains: a data.frame containing MCMC outputs (which are also
#' stored in the file indicated in \code{res.file.name}).
#'
#' \item collec.dates: (data) the collection dates.
#'
#' \item w: (data) the generation time distribution (argument \code{w.dens})
#'
#' \item f: (data) the distribution of the time to collection (argument
#' \code{f.dens})
#'
#' \item D: a matrix of genetic distances (in number of mutations) between all
#' pairs of sequences.
#'
#' \item idx.dna: (data) the index of the case each dna sequence corresponds to
#'
#' \item tune.end: an integer indicating at which iteration the proposal
#' auto-tuning procedures all stopped.
#'
#' \item find.import: a logical indicating if imported cases were to be
#' automatically detected.
#'
#' \item burnin: an integer indicating the pre-defined burnin, used when
#' detecting imported cases.
#'
#' \item find.import.at: an integer indicating at which iteration of the
#' preliminary MCMC imported cases were detected.
#'
#' \item n.runs: the number of independent runs used.
#'
#' \item call: the matched call.  }
#' @author Thibaut Jombart (\email{t.jombart@@imperial.ac.uk})
#' @seealso \itemize{ \item \link{plotChains} to visualize MCMC chains.
#'
#' \item \link{transGraph} and \link{get.tTree} to represent transmission
#' trees.
#'
#' \item \link{get.R} and \link{get.Rt} to get reproduction numbers
#' distributions.
#'
#' \item \link{get.incid} to get estimates of incidence.
#'
#' \item \link{get.mu} to get the mutation rate distribution.
#'
#' \item \link{simOutbreak} to simulate outbreaks.
#'
#' \item \link{selectChains} to select chains from parallel runs which
#' converged towards different posterior modes.
#'
#' \item \link{fakeOutbreak}, a toy dataset used to illustrate the method.
#'
#' \item For more resources including tutorials, forums, etc., see:
#' \url{http://sites.google.com/site/therepiproject/r-pac/outbreaker}
#'
#' }
#' @references Jombart T, Cori A, Didelot X, Cauchemez S, Fraser C and Ferguson
#' N (accepted).  Bayesian reconstruction of disease outbreaks by combining
#' epidemiologic and genomic data. PLoS Computational Biology.
#' @examples
#'
#'
#' ## EXAMPLE USING TOYOUTBREAK ##
#' ## LOAD DATA, SET RANDOM SEED
#' data(fakeOutbreak)
#' attach(fakeOutbreak)
#'
#' ## VISUALIZE DYNAMICS
#' matplot(dat$dynam, type="o", pch=20, lty=1,
#'    main="Outbreak dynamics", xlim=c(0,28))
#' legend("topright", legend=c("S","I","R"), lty=1, col=1:3)
#'
#' ## VISUALIZE TRANSMISSION TREE
#' plot(dat, annot="dist", main="Data - transmission tree")
#' mtext(side=3, "arrow annotations are numbers of mutations")
#'
#'
#' \dontrun{
#' ## RUN OUTBREAKER - PARALLEL VERSION
#' ## (takes < 1 min))
#' set.seed(1)
#' res <-  outbreaker.parallel(n.runs=4, dna=dat$dna,
#'    dates=collecDates,w.dens=w, n.iter=5e4)
#' }
#'
#'
#' ## ASSESS CONVERGENCE OF CHAINS
#' plotChains(res)
#' plotChains(res, burnin=2e4)
#'
#' ## REPRESENT POSTERIOR ANCESTRIES
#' transGraph(res, annot="", main="Posterior ancestries", thres=.01)
#'
#' ## GET CONSENSUS ANCESTRIES
#' tre <- get.tTree(res)
#' plot(tre, annot="", main="Consensus ancestries")
#'
#' ## SHOW DISCREPANCIES
#' col <- rep("lightgrey", 30)
#' col[which(dat$ances != tre$ances)] <- "pink"
#' plot(tre, annot="", vertex.color=col, main="Consensus ancestries")
#' mtext(side=3, text="cases with erroneous ancestries in pink")
#'
#' ## GET EFFECTIVE REPRODUCTION OVER TIME
#' get.Rt(res)
#'
#' ## GET INDIVIDUAL EFFECTIVE REPRODUCTION
#' head(get.R(res))
#' boxplot(get.R(res), col="grey", xlab="Case",
#'         ylab="Effective reproduction number")
#'
#' ## GET MUTATION RATE PER TIME UNIT
#' ## per genome
#' head(get.mu(res))
#'
#' ## per nucleotide
#' mu <- get.mu(res, genome.size=1e4)
#' head(mu)
#'
#' summary(mu)
#' hist(mu, border="lightgrey", col="grey", xlab="Mutation per day and nucleotide",
#'      main="Posterior distribution of mutation rate")
#'
#' detach(fakeOutbreak)
#'
#'
#' @importFrom stats dexp
#'
outbreaker <- function(dates, dna=NULL,
                       w.dens, f.dens=w.dens,
                       init.tree=c("seqTrack","star","random"),
                       init.mu=NULL,
                       move.ances=TRUE, move.Tinf=TRUE, move.mut=TRUE,
                       n.iter=1e5, sample.every=500, tune.every=500){

    ## CHECKS / PROCESS DATA ##

    ## DATES ##
    ## conversions
    if(inherits(dates, "Date")) dates <- dates-min(dates)
    if(inherits(dates, "POSIXct")) dates <- difftime(dates, min(dates), units="days")
    dates <- as.integer(round(dates))

    ## global variable (nb of cases)
    N <- length(dates)

    ## DNA SEQUENCES ##
    ## check type of input ##
    if(!inherits(dna, "DNAbin")) stop("dna is not a DNAbin object.")
    if(!is.matrix(dna)) dna <- as.matrix(dna)
    if(is.character(dates)) stop("dates are characters; they must be integers or dates with Date format (see ?as.Date)")

    ## global variables
    L <- ncol(dna) #  (genome length)
    D <- as.matrix(dist.dna(dna, model="N")) # distance matrix


    ## DENSITIES ##
    ## checks
    if(any(w.dens)<0) {
        stop("w.dens has negative entries (these should be probabilities!)")
    }
    if(any(f.dens)<0) {
        stop("f.dens has negative entries (these should be probabilities!)")
    }

    ## set p(T=0) to zero
    w.dens[1] <- f.dens[1] <- 0

    ## standardize densities
    w.dens <- w.dens/sum(w.dens)
    f.dens <- f.dens/sum(f.dens)

    ## find range
    max.range <- diff(range(dates))

    ## add an exponential tail summing to 1e-4 to 'w'
    ## to cover the span of the outbreak
    ## (avoids starting with -Inf temporal loglike)
    if(length(w.dens)<max.range) {
        length.to.add <- (max.range-length(w.dens)) + 10 # +10 to be on the safe side
        val.to.add <- dexp(1:length.to.add, 1)
        val.to.add <- 1e-4*(val.to.add/sum(val.to.add))
        w.dens <- c(w.dens, val.to.add)
        w.dens <- w.dens/sum(w.dens)
    }

    ## get log-densities
    log.w.dens <- log(w.dens)
    log.f.dens <- log(f.dens)


    ## INITIAL PARAMETER VALUES ##
    ## TREE ##
    if(is.character(init.tree)) {
        init.tree <- match.arg(init.tree)
    } else {
        if(length(init.tree) != length(dates)) stop("inconvenient length for init.tree")
        init.tree[init.tree<1 | init.tree>N] <- NA
    }

    ## MUTATION RATE 'MU' ##
    if(!is.null(init.mu) && init.mu1<0) stop("init.mu1 < 0")


    ## complete w.dens ##
    max.range <- diff(range(dates))
    ## add an exponential tail summing to 1e-4 to 'w'
    ## to cover the span of the outbreak
    ## (avoids starting with -Inf temporal loglike)
    if(length(w.dens)<max.range) {
        length.to.add <- (max.range-length(w.dens)) + 10 # +10 to be on the safe side
        val.to.add <- dexp(1:length.to.add, 1)
        val.to.add <- 1e-4*(val.to.add/sum(val.to.add))
        w.dens <- c(w.dens, val.to.add)
        w.dens <- w.dens/sum(w.dens)
    }



    ## find initial tree ##
    ## get temporal ordering constraint:
    ## canBeAnces[i,j] is 'i' can be ancestor of 'j'
    canBeAnces <- outer(dates,dates,FUN="<") # strict < is needed as we impose w(0)=0
    diag(canBeAnces) <- FALSE

    if(is.character(init.tree)){
        ## check init
        if(init.tree=="seqTrack" && !all(1:n.ind %in% idx.dna)) {
            warning("Can't use seqTrack initialization with missing DNA sequences - using a star-like tree")
            init.tree <- "star"
        }

        ## seqTrack init
        if(init.tree=="seqTrack"){
            D.temp <- D
            D.temp[!canBeAnces] <- 1e30
            ances <- apply(D.temp,2,which.min)
            ances[dates==min(dates)] <- NA
            ances <- as.integer(ances)
        }

        ## star-shape init
        if(init.tree=="star"){
            ances <- rep(which.min(dates), length(dates))
            ances[dates==min(dates)] <- 0
            ances <- as.integer(ances-1) # put on C scale
        }

        ## random init
        if(init.tree=="random"){
            ances <- apply(canBeAnces, 2, function(e) ifelse(length(which(e))>0, sample(which(e),1), NA) )
            ances <- ances-1
            ances[is.na(ances)] <- -1L
            ances <- as.integer(ances)
        }
    }


    ## MCMC ##
    ## create output templates ##
    out.post <- out.prior <- out.ll <- double(n.iter)
    mu <- init.mu


    ## initialize algorithm ##
    out.ll[1] <- ll.all(times=dates, ances=ances, log.w=log.w.dens, D=D, mu=mu, gen.length=L)

    out <- NULL
    return(out)
} # end outbreaker







## #' @rdname outbreaker
## #' @export
## outbreaker.parallel <- function(n.runs, parallel=TRUE, n.cores=NULL,
##                                 dna=NULL, dates, idx.dna=NULL, mut.model=1, spa.model=0,
##                                 w.dens, f.dens=w.dens,
##                                 dist.mat=NULL,
##                                 init.tree=c("seqTrack","random","star"),
##                                 init.kappa=NULL,
##                                 init.mu1=NULL, init.mu2=init.mu1, init.spa1=NULL,
##                                 n.iter=1e5, sample.every=500, tune.every=500,
##                                 burnin=2e4, import.method=c("genetic","full","none"),
##                                 find.import.n=50,
##                                 pi.prior1=10, pi.prior2=1, spa1.prior=1,
##                                 move.mut=TRUE, move.ances=TRUE, move.kappa=TRUE,
##                                 move.Tinf=TRUE, move.pi=TRUE, move.spa=TRUE,
##                                 outlier.threshold = 5, max.kappa=10,
##                                 quiet=TRUE, res.file.name="chains.txt", tune.file.name="tuning.txt", seed=NULL){

##     ## SOME CHECKS ##
##     if(parallel && is.null(n.cores)){
##         n.cores <- detectCores()
##         n.cores <- min(n.cores, 6)
##     }


##     ## GET FILE NAMES ##
##     res.file.names <- paste("run", 1:n.runs, "-", res.file.name, sep="")
##     tune.file.names <- paste("run", 1:n.runs, "-", tune.file.name, sep="")


##     ## HANDLE SEED ##
##     if(is.null(seed)){
##         seed <- as.integer(runif(n.runs,min=0,max=2e9))
##     } else {
##         seed <- rep(seed, length=n.runs)
##     }


##     ## COMPUTATIONS ##
##     if(parallel){
##         ## create cluster ##
##         clust <- makeCluster(n.cores)

##         ## load outbreaker for each child ##
##         clusterEvalQ(clust, library(outbreaker))

##         ## transfer data onto each child ##
##         listArgs <- c("dna", "dates", "idx.dna", "mut.model", "spa.model", "w.dens", "f.dens", "dist.mat", "init.tree", "init.kappa", "n.iter",
##                       "sample.every", "tune.every", "burnin", "import.method", "find.import.n", "pi.prior1", "pi.prior2", "init.mu1", "init.mu2",
##                       "init.spa1", "move.mut", "spa1.prior", "move.mut", "move.ances", "move.kappa", "move.Tinf", "move.pi", "move.spa",
##                       "outlier.threshold", "max.kappa", "res.file.names", "tune.file.names", "seed")

##         clusterExport(clust, listArgs, envir=environment())

##         ## set calls to outbreaker on each child ##
##         res <- parLapply(clust, 1:n.runs, function(i)  outbreaker(dna=dna, dates=dates, idx.dna=idx.dna,
##                                                                   mut.model=mut.model, spa.model=spa.model,
##                                                                   w.dens=w.dens,
##                                                                   f.dens=f.dens,
##                                                                   dist.mat=dist.mat, ## locations=locations,
##                                                                   init.tree=init.tree, init.kappa=init.kappa,
##                                                                   n.iter=n.iter, sample.every=sample.every,
##                                                                   tune.every=tune.every, burnin=burnin,
##                                                                   import.method=import.method,
##                                                                   find.import.n=find.import.n,
##                                                                   pi.prior1=pi.prior1, pi.prior2=pi.prior2,
##                                                                   spa1.prior=spa1.prior,
##                                                                   init.mu1=init.mu1, init.mu2=init.mu2, init.spa1=init.spa1,
##                                                                   move.mut=move.mut, move.ances=move.ances, move.kappa=move.kappa,
##                                                                   move.Tinf=move.Tinf, move.pi=move.pi, move.spa=move.spa,
##                                                                   outlier.threshold = outlier.threshold, max.kappa=max.kappa,
##                                                                   quiet=TRUE, res.file.name=res.file.names[i],
##                                                                   tune.file.name=tune.file.names[i], seed=seed[i]))

##         ## close parallel processes ##
##         stopCluster(clust)

##         ## Version with mclapply - doesn't work on windows ##
##         ## res <- mclapply(1:n.runs, function(i)  outbreaker(dna=dna, dates=dates, idx.dna=idx.dna, w.dens=w.dens, w.trunc=w.trunc,
##         ##                                                     init.tree=init.tree, init.kappa=init.kappa,
##         ##                                                     n.iter=n.iter, sample.every=sample.every,
##         ##                                                     tune.every=tune.every, burnin=burnin,
##         ##                                                     find.import=find.import, find.import.n=find.import.n,
##         ##                                                     pi.prior1=pi.prior1, pi.prior2=pi.prior2,
##         ##                                                     init.mu1=init.mu1, init.mu2=init.mu2,
##         ##                                                     move.mut=move.mut, move.ances=move.ances, move.kappa=move.kappa,
##         ##                                                     move.Tinf=move.Tinf, move.pi=move.pi,
##         ##                                                     quiet=TRUE, res.file.name=res.file.names[i],
##         ##                                                     tune.file.name=tune.file.names[i], seed=seed[i]),
##         ##                   mc.cores=n.cores, mc.silent=FALSE, mc.cleanup=TRUE, mc.preschedule=TRUE, mc.set.seed=TRUE)
##     } else {
##         res <- lapply(1:n.runs, function(i)  outbreaker(dna=dna, dates=dates, idx.dna=idx.dna,
##                                                         mut.model=mut.model, spa.model=spa.model,
##                                                         w.dens=w.dens,
##                                                         f.dens=f.dens,
##                                                         dist.mat=dist.mat,
##                                                         init.tree=init.tree, init.kappa=init.kappa,
##                                                         n.iter=n.iter, sample.every=sample.every,
##                                                         tune.every=tune.every, burnin=burnin,
##                                                         import.method=import.method,
##                                                         find.import.n=find.import.n,
##                                                         pi.prior1=pi.prior1, pi.prior2=pi.prior2,
##                                                         spa1.prior=spa1.prior,
##                                                         init.mu1=init.mu1, init.mu2=init.mu2, init.spa1=init.spa1,
##                                                         move.mut=move.mut, move.ances=move.ances, move.kappa=move.kappa,
##                                                         move.Tinf=move.Tinf, move.pi=move.pi, move.spa=move.spa,
##                                                         outlier.threshold = outlier.threshold, max.kappa=max.kappa,
##                                                         quiet=TRUE, res.file.name=res.file.names[i],
##                                                         tune.file.name=tune.file.names[i], seed=seed[i]))
##     }


##     ## MERGE RESULTS ##
##     res.old <- res
##     res <- res[[1]]
##     res$tune.end <- max(sapply(res.old, function(e) e$tune.end))
##     res$chains <- Reduce(rbind, lapply(res.old, function(e) e$chains))
##     res$chains$run <- factor(rep(1:n.runs, each=nrow(res.old[[1]]$chains)))
##     res$n.runs <- n.runs
##     res$call <- match.call()

##     ## RETURN RESULTS ##
##     return(res)
## } # end outbreaker.parallel

