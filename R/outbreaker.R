
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
#' @examples
#' \dontrun{
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
#' }
#'
#' @seealso outbreaker.config to see default parameters / set parameters
#'
#' @importFrom stats dexp rnorm runif
#' @importFrom coda mcmc
#'
outbreaker <- function(dates, dna=NULL,
                       w.dens, f.dens=w.dens,
                       config=outbreaker.config()){

    ## CHECKS / PROCESS DATA ##
    data <- check.data(dates, dna, w.dens, f.dens)







    ## MCMC ##
    ## create output templates ##
    out.size <- round(n.iter/sample.every) + 1 # +1 because initial state is always there
    out.post <- out.prior <- out.like <- out.mu <- double(out.size)
    out.ances <- as.list(integer(out.size))
    out.t.inf <- as.list(integer(out.size))

    ## initialize algorithm and outputs ##
    current.mu <- out.mu[1] <- init.mu
    current.ances <- out.ances[[1]] <- ances
    current.t.inf <- out.t.inf[[1]] <- dates - which.max(f.dens) + 1
    out.like[1] <- ll.all(t.inf=out.t.inf[[1]], ances=out.ances[[1]],
                          log.w=log.w.dens, log.f=log.f.dens,
                          sampling.times=dates,
                          D=D, mu=init.mu, gen.length=L)
    out.prior[1] <- prior.all(init.mu)
    out.post[1] <- out.like[1] + out.prior[1]
    OUT.COUNTER <- 1

    ## initialize pre-drawn random arrays ##
    RAND.MU <- rnorm(n.iter, mean=0, sd=sd.mu)
    RAND.ACC.MU <- log(runif(n.iter))
    RAND.ACC.T.INF <- log(runif(n.iter))
    RAND.ACC.ANCES <- log(runif(n.iter))

    ## run MCMC ##
    for(i in 2:n.iter){
        ## move infection dates ##
        if(move.t.inf){
            current.t.inf <- move.t.inf(sampling.times=dates, log.w=log.w.dens, log.f=log.f.dens,
                                        t.inf=current.t.inf, ances=current.ances,
                                        lunif=RAND.ACC.T.INF[i])
        }

        ## move ancestries ##
        if(move.ances){
            current.ances <- move.ances(sampling.times=dates, D=D,
                                        gen.length=L, log.w=log.w.dens, log.f=log.f.dens,
                                        t.inf=current.t.inf, ances=current.ances, mu=current.mu,
                                        lunif=RAND.ACC.ANCES[i])
        }

        ## move mu ##
        if(move.mu){
            current.mu <- move.mu(D=D, gen.length=L, ances=current.ances, mu=current.mu,
                                  dev=RAND.MU[i], lunif=RAND.ACC.MU[i])
        }

        ## store outputs if needed
        if((i %% sample.every) == 0){
            OUT.COUNTER <- OUT.COUNTER + 1
            ## store like, prior, post
            out.like[OUT.COUNTER] <- ll.all(t.inf=current.t.inf, ances=current.ances,
                                            log.w=log.w.dens, log.f=log.f.dens,
                                            sampling.times=dates,
                                            D=D, mu=init.mu, gen.length=L)
            out.prior[OUT.COUNTER] <- prior.all(current.mu)
            out.post[OUT.COUNTER] <- out.like[OUT.COUNTER] + out.prior[OUT.COUNTER]

            ## parameters and augmented data
            out.mu[OUT.COUNTER] <- current.mu
            out.ances[[OUT.COUNTER]] <- current.ances
            out.t.inf[[OUT.COUNTER]] <- current.t.inf
        }

    } # end of the chain


    ## SHAPE RESULTS ##
    out.ances <- matrix(unlist(out.ances), ncol=N, byrow=TRUE)
    colnames(out.ances) <- paste("alpha",1:N, sep=".")
    out.t.inf <- matrix(unlist(out.t.inf), ncol=N, byrow=TRUE)
    colnames(out.t.inf) <- paste("t.inf",1:N, sep=".")
    out <- data.frame(post=out.post, like=out.like, prior=out.prior, mu=out.mu, out.ances, out.t.inf)
    out <- coda::mcmc(out)

    ## RETURN ##
    return(out)
} # end outbreaker

