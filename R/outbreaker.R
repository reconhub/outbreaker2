
#' Outbreaker2: disease outbreak reconstruction using epidemiological and genetic data
#'
#'
#' @export
#'
#' @aliases outbreaker
#'
#' @rdname outbreaker
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
#' @param init.tree the tree used to initialize the MCMC. Can be either a
#' character string indicating how this tree should be computed, or a vector of
#' integers corresponding to the tree itself, where the i-th value corresponds
#' to the index of the ancestor of 'i' (i.e., \code{init.tree[i]} is the
#' ancestor of case \code{i}). Accepted character strings are "seqTrack" (uses
#' seqTrack output as initialize tree), "random" (ancestor randomly selected
#' from preceding cases), and "star" (all cases coalesce to the first case).
#' Note that for SeqTrack, all cases should have been sequenced.
#'
#' @param n.iter an integer indicating the number of iterations in the MCMC,
#' including the burnin period; defaults to \code{100,000}.
#'
#' @param init.mu initial values for the mutation rates
#'
#' @param move.mut logical indicating whether the mutation rates
#' should be estimated ('moved' in the MCMC), or not, all defaulting to TRUE.
#'
#' @param move.ances,move.kappa,move.Tinf vectors of logicals of length 'n'
#' indicating for which cases different components should be moved during the
#' MCMC.
#'
#' @author Thibaut Jombart (\email{t.jombart@@imperial.ac.uk})
#'
#' @references Jombart T, Cori A, Didelot X, Cauchemez S, Fraser C and Ferguson
#' N (2014).  Bayesian reconstruction of disease outbreaks by combining
#' epidemiologic and genomic data. PLoS Computational Biology.
#'
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
#' @importFrom stats dexp rnorm runif
#' @importFrom coda mcmc
#'
outbreaker <- function(dates, dna=NULL,
                       w.dens, f.dens=w.dens,
                       init.tree=c("seqTrack","star","random"),
                       init.mu=1e-4,
                       move.ances=TRUE, move.Tinf=TRUE, move.mut=TRUE,
                       n.iter=10, sd.mu=0.0001){

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
    if(any(w.dens<0)) {
        stop("w.dens has negative entries (these should be probabilities!)")
    }
    if(any(f.dens<0)) {
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
    if(!is.null(init.mu) && init.mu<0) stop("init.mu < 0")


    ## FIND INITIAL TREE ##
    ## get temporal ordering constraint:
    ## canBeAnces[i,j] is 'i' can be ancestor of 'j'
    canBeAnces <- outer(dates,dates,FUN="<") # strict < is needed as we impose w(0)=0
    diag(canBeAnces) <- FALSE

    if(is.character(init.tree)){
        ## check init
        if(init.tree=="seqTrack" && is.null(dna)) {
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
        }

        ## random init
        if(init.tree=="random"){
            ances <- rances(dates)
        }
    }


    ## MCMC ##
    ## create output templates ##
    out.post <- out.prior <- out.like <- out.mu <- double(n.iter)
    out.ances <- as.list(integer(n.iter))
    out.t.inf <- as.list(integer(n.iter))

    ## initialize algorithm and outputs ##
    out.mu[1] <- init.mu
    out.ances[[1]] <- ances
    out.t.inf[[1]] <- dates - which.max(f.dens) + 1
    out.like[1] <- ll.all(t.inf=out.t.inf[[1]], ances=out.ances[[1]],
                          log.w=log.w.dens, log.f=log.f.dens,
                          sampling.times=dates,
                          D=D, mu=init.mu, gen.length=L)
    out.prior[1] <- prior.all(init.mu)
    out.post[1] <- out.like[1] + out.prior[1]

    ## initialize pre-drawn random arrays ##
    RAND.MU <- rnorm(n.iter, mean=0, sd=sd.mu)
    RAND.ACC.MU <- log(runif(n.iter))
    RAND.ACC.T.INF <- log(runif(n.iter))
    RAND.ACC.ANCES <- log(runif(n.iter))

    ## run MCMC ##
    for(i in 2:n.iter){
        ## move ancestries ##
        out.ances[[i]] <- move.ances(t.inf=t.inf, sampling.times=sampling.times, D=D,
                                     gen.length=gen.length, log.w=log.w, log.f=log.f,
                                     ances=out.ances[[i-1]], mu=mu, lunif=RAND.ACC.ANCES)

        ## move infection dates ##
        out.t.inf[[i]] <- move.t.inf(t.inf=out.t.inf[[i-1]], log.w=log.w.dens, log.f=log.f.dens,
                                     ances=out.ances[[i]], sampling.times=dates,
                                     lunif=RAND.ACC.T.INF[i])

        ## move mu ##
        out.mu[i] <- move.mu(D=D, gen.length=L, ances=out.ances[[i]], mu=out.mu[i-1],
                             dev=RAND.MU[i], lunif=RAND.ACC.MU[i])

    } # end of the chain


    ## SHAPE RESULTS ##
    out <- data.frame(post=out.post, like=out.like, prior=out.prior, mu=out.mu)
    out <- coda::mcmc(out)

    ## RETURN ##
    return(out)
} # end outbreaker

