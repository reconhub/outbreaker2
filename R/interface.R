
##################
## main functions
##################
outbreaker <- function(dna=NULL, dates, idx.dna=NULL,
                       mut.model=1, spa.model=0,
                       w.dens, w.trunc=length(w.dens),
                       f.dens=w.dens, f.trunc=length(f.dens),
                       dist.mat=NULL,
                       init.tree=c("seqTrack","random","star"),
                       init.kappa=NULL, init.mu1=NULL, init.mu2=init.mu1, init.spa1=NULL,
                       n.iter=1e5, sample.every=500, tune.every=500,
                       burnin=2e4, import.method=c("genetic","full","none"),
                       find.import.n=50,
                       pi.prior1=10, pi.prior2=1, spa1.prior=1,
                       move.mut=TRUE, move.ances=TRUE, move.kappa=TRUE,
                       move.Tinf=TRUE, move.pi=TRUE, move.spa=TRUE,
                       outlier.threshold = 5, max.kappa=10,
                       quiet=TRUE, res.file.name="chains.txt",
                       tune.file.name="tuning.txt", seed=NULL){

    ## CHECKS ##
    ## if(!require(ape)) stop("the ape package is required but not installed")
    ## RE-ORDERING OF IMPORT METHOD TO MATCH C SIDE:
    ## none:0L
    ## genetic: 1L
    ## full: 2L
    import.method <- match.arg(import.method)
    import.method <- as.integer(match(import.method, c("none", "genetic","full")))-1L


    ## HANDLE MISSING DNA ##
    useDna <- !is.null(dna)
    if(is.null(dna)){
        dna <- as.DNAbin(matrix('a',ncol=10,nrow=length(dates)))
        move.mut <- FALSE
        mut.model <- 0L
        if(is.character(init.tree) && match.arg(init.tree)=="seqTrack") init.tree <- "star"
        init.mu1 <- init.mu2 <-0
        init.gamma <- 1
    }

    if(!inherits(dna, "DNAbin")) stop("dna is not a DNAbin object.")
    if(!is.matrix(dna)) dna <- as.matrix(dna)
    if(is.character(dates)) stop("dates are characters; they must be integers or dates with POSIXct format (see ?as.POSIXct)")
    if(is.character(init.tree)) {
        init.tree <- match.arg(init.tree)
    } else {
        if(length(init.tree) != length(dates)) stop("inconvenient length for init.tree")
        init.tree[is.na(init.tree)|init.tree<1] <- 0
        if(max(init.tree)>length(dates)) stop("inconvenient values in init.tree (some indices > n)")
        ances <- as.integer(init.tree-1) # translate indices on C scale (0:(n-1))
    }
    if(length(w.dens)<w.trunc) stop(paste("incomplete w.dens: values needed from t=0 to t=", w.trunc-1,sep=""))
    w.dens[1] <- 0 # force w_0 = 0
    w.dens[w.dens<0] <- 0
    if(sum(w.dens) <= 1e-14) stop("w.dens is zero everywhere")
    if(!is.null(init.mu1) && init.mu1<0) stop("init.mu1 < 0")
    if(!is.null(init.mu2) && init.mu2<0) stop("init.mu2 < 0")


    ## PROCESS INPUTS ##
    ## dna ##
    n.seq <- as.integer(nrow(dna))
    n.ind <- as.integer(length(dates))
    n.nucl <- as.integer(ncol(dna))
    dnaraw <- unlist(as.list(dna),use.names=FALSE)
    ## if(n.ind != length(dates)) stop(paste("dna and dates have different number of individuals -",n.ind,"versus",length(dates)))

    ## handle dates ##
    if(is.numeric(dates)){
        if(sum(abs(dates-round(dates))>1e-15)) warning("dates have been rounded to nearest integers")
        dates <- as.integer(round(dates))
    }

    if(inherits(dates, "POSIXct")){
        dates <- difftime(dates, min(dates), units="days")
    }
    dates <- as.integer(dates)


    ## handle idx.dna ##
    ## need to go from: id of case for each sequence (idx.dna)
    ## to: position of the sequence in DNA matrix for each case
    ## -1 is used for missing sequences
    if(is.null(idx.dna)) {
        idx.dna <- 1:nrow(dna)
    }

    if(any(!idx.dna %in% 1:n.ind)) stop("DNA sequences provided for unknown cases (some idx.dna not in 1:n.ind)")
    if(length(idx.dna)!=nrow(dna)) stop("length of idx.dna does not match the number of DNA sequences")

    idx.dna.for.cases <- match(1:n.ind, idx.dna)
    idx.dna.for.cases[is.na(idx.dna.for.cases)] <- 0
    idx.dna.for.cases <- as.integer(idx.dna.for.cases-1) # for C

    ## check mutational model ##
    ## check model
    mut.model <- as.integer(mut.model)
    if(!mut.model %in% c(0L,1L,2L)) stop("unknown mutational model requested; accepted values are: 0, 1, 2")
    ## model 0: no evolution model
    if(mut.model==0L){
        init.gamma <- 1
        init.mu1 <- 1
        move.mut <- FALSE
    }
    ## determine gamma
    if(!is.null(init.mu1) && !is.null(init.mu2)){
        init.gamma <- init.mu2/init.mu1
        if(is.na(init.gamma) || is.infinite(init.gamma)) init.gamma <- 1 # in case rates are both 0
    } else{
        init.gamma <- 1
    }
    ## force gamma to 1 for model 1
    if(mut.model==1L){
        init.gamma <- 1
    }

    ## check generation time function ##
    w.dens <- as.double(w.dens)
    w.dens <- w.dens/sum(w.dens)
    if(any(is.na(w.dens))) stop("NAs in w.dens after normalization")
    w.trunc <- as.integer(w.trunc)

    ## check collection time function ##
    f.dens <- as.double(f.dens)
    f.dens <- f.dens/sum(f.dens)
    if(any(is.na(f.dens))) stop("NAs in f.dens after normalization")
    f.trunc <- as.integer(f.trunc)

    ## check spatial distances ##
    if(!is.null(dist.mat)){
        if(!inherits(dist.mat,"matrix")) dist.mat <- as.matrix(dist.mat)
        if(nrow(dist.mat) != ncol(dist.mat)) stop("matrix of distances (dist.mat) is not square")
        if(nrow(dist.mat) != length(dates)) stop("wrong dimension for the matrix of distances")
        if(any(is.na(dist.mat))) stop("NAs in the distance matrix")
    } else {
        if(spa.model>0) stop("spatial model requested but dist.mat not provided")
        dist.mat <- matrix(0, ncol=length(dates), nrow=length(dates))
        spa.model <- 0L
    }

    ## check spatial model ##
    spa.model <- as.integer(spa.model)
    if(!spa.model %in% c(0L, 1L, 2L)) stop("unknown spatial model requested; accepted values are: 0, 1, 2")
    ## model 0: no spatial info
    if(spa.model == 0L) {
        dist.mat <- matrix(0, ncol=length(dates), nrow=length(dates))
        init.spa1 <- init.spa2 <- 0
    }
    ## model 1: normal dispersal
    if(spa.model > 0L) {
        if(is.null(init.spa1)) init.spa1 <- 1
        ## if(is.null(init.spa2)) init.spa2 <- 0
        init.spa2 <- 0
        spa1.prior <- max(0.0, spa1.prior)
    }
    ## model 2: stratified dispersal
    if(spa.model == 2L){
        stop("Only available spatial models are 0 (none) and 1 (exponential diffusion)")
        ## if(is.null(locations)) stop("Spatial model 2 needs locations for stratified dispersal")
        ## if(length(locations)!=length(dates)) stop("wrong length for argument 'locations'")
        ## if(is.character(locations)) locations <- factor(locations)
        ## locations <- as.integer(locations)
    }


    ## init.kappa ##
    ## if NULL, will be ML assigned (code is kappa_i<0)
    if(is.null(init.kappa)) init.kappa <- rep(0L,n.ind)
    init.kappa <- as.integer(rep(init.kappa, length=n.ind)) #recycle


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
            D <- as.matrix(dist.dna(dna, model="TN93"))
            D[!canBeAnces] <- 1e15
            ances <- apply(D,2,which.min)-1 # -1 for compatibility with C
            ances[dates==min(dates)] <- -1 # unknown ancestor
            ances <- as.integer(ances)
        }

        ## random init
        if(init.tree=="random"){
            ances <- apply(canBeAnces, 2, function(e) ifelse(length(which(e))>0, sample(which(e),1), NA) )
            ances <- ances-1
            ances[is.na(ances)] <- -1L
            ances <- as.integer(ances)
        }

        ## star-shape init
        if(init.tree=="star"){
            ances <- rep(which.min(dates), length(dates))
            ances[dates==min(dates)] <- 0
            ances <- as.integer(ances-1) # put on C scale
        }

        ## ## no ancestry init
        ## if(init.tree=="none"){
        ##     ances <- as.integer(rep(-1,length(dates)))
        ## }
    }

    ## handle seed ##
    if(is.null(seed)){
        seed <- as.integer(runif(1,min=0,max=2e9))
    }

    ## handle import.method ##
    if(mut.model==0L && import.method==1L) import.method <- 2L

    ## handle find.import ##
    find.import <- import.method > 0L
    if(find.import){
        find.import.n <- max(find.import.n,30) # import at least based on 30 values
        find.import.at <- as.integer(round(burnin + find.import.n*sample.every))
        if(find.import.at>=n.iter) stop(paste("n.iter (", n.iter, ") is less than find.import.at (", find.import.at,")", sep=""))
    } else {
        find.import.at <- as.integer(0)
    }


    ## coerce type for remaining arguments ##
    n.iter <- as.integer(n.iter)
    sample.every <- as.integer(sample.every)
    tune.every <- as.integer(tune.every)
    pi.prior1 <- as.double(pi.prior1)
    pi.prior2 <- as.double(pi.prior2)
    ## phi.param1 <- as.double(phi.param1)
    ## phi.param2 <- as.double(phi.param2)
    phi.param1 <- phi.param2 <- as.double(1)
    if(is.null(init.mu1)) {
        init.mu1 <- 0.5/ncol(dna)
    }
    init.mu1 <- as.double(init.mu1)
    init.gamma <- as.double(init.gamma)
    init.spa1 <- as.double(init.spa1)
    ## init.spa2 <- as.double(init.spa2)
    init.spa2 <- as.double(1)
    spa1.prior <- as.double(spa1.prior)
    ## spa2.prior <- as.double(spa2.prior)
    spa2.prior <- as.double(1)
    move.mut <- as.integer(move.mut)
    move.ances <- as.integer(rep(move.ances, length=n.ind))
    move.kappa <- as.integer(rep(move.kappa, length=n.ind))
    move.Tinf <- as.integer(move.Tinf)
    move.pi <- as.integer(move.pi)
    ## move.phi <- as.integer(move.phi)
    move.phi <- 0L
    move.spa <- as.integer(move.spa)
    quiet <- as.integer(quiet)
    res.file.name <- as.character(res.file.name)[1]
    tune.file.name <- as.character(tune.file.name)[1]
    burnin <- as.integer(burnin)
    outlier.threshold <- as.double(outlier.threshold)
    max.kappa <- as.integer(max.kappa)
    ##locations <- as.integer(locations)
    locations <- rep(0L, length(dates))


    ## create empty output vector for genetic distances ##
    dna.dist <- integer(n.ind*(n.ind-1)/2)
    stopTuneAt <- integer(1)

    temp <- .C("R_outbreaker",
               dnaraw, dates, n.ind, n.seq, n.nucl,  idx.dna.for.cases, mut.model,
               w.dens, w.trunc, f.dens, f.trunc,
               dist.mat, locations, spa.model,
               ances, init.kappa, n.iter, sample.every, tune.every,
               pi.prior1, pi.prior2, phi.param1, phi.param2, init.mu1, init.gamma,
               init.spa1, init.spa2, spa1.prior, spa2.prior,
               move.mut, move.ances, move.kappa, move.Tinf,
               move.pi, move.phi, move.spa,
               import.method, find.import.at, burnin, outlier.threshold,
               max.kappa, quiet,
               dna.dist, stopTuneAt, res.file.name, tune.file.name, seed,
               PACKAGE="outbreaker")

    D <- temp[[43]]
    D[D<0] <- NA
    stopTuneAt <- temp[[44]]

    cat("\nComputations finished.\n\n")

    ## make D a 'dist' object ##
    attr(D,"Size") <- n.ind
    attr(D,"Diag") <- FALSE
    attr(D,"Upper") <- FALSE
    class(D) <- "dist"


    ## BUILD OUTPUT ##
    ## read table
    chains <- read.table(res.file.name, header=TRUE, stringsAsFactors=FALSE,
                         colClasses=c("integer", rep("numeric",7+n.ind*2)))

    chains$run <- rep(1, nrow(chains))
    call <- match.call()
    res <- list(chains=chains, collec.dates=dates, w=w.dens[1:w.trunc], f=f.dens[1:f.trunc], D=D, idx.dna=idx.dna, tune.end=stopTuneAt,
                burnin=burnin, import.method=import.method, find.import.at=find.import.at, n.runs=1, call=call)

    return(res)
} # end outbreaker









###############################
## version with multiple runs
###############################
outbreaker.parallel <- function(n.runs, parallel=TRUE, n.cores=NULL,
                                dna=NULL, dates, idx.dna=NULL, mut.model=1, spa.model=0,
                                w.dens, w.trunc=length(w.dens),
                                f.dens=w.dens, f.trunc=length(f.dens),
                                dist.mat=NULL,
                                init.tree=c("seqTrack","random","star"),
                                init.kappa=NULL,
                                init.mu1=NULL, init.mu2=init.mu1, init.spa1=NULL,
                                n.iter=1e5, sample.every=500, tune.every=500,
                                burnin=2e4, import.method=c("genetic","full","none"),
                                find.import.n=50,
                                pi.prior1=10, pi.prior2=1, spa1.prior=1,
                                move.mut=TRUE, move.ances=TRUE, move.kappa=TRUE,
                                move.Tinf=TRUE, move.pi=TRUE, move.spa=TRUE,
                                outlier.threshold = 5, max.kappa=10,
                                quiet=TRUE, res.file.name="chains.txt", tune.file.name="tuning.txt", seed=NULL){

    ## SOME CHECKS ##
    if(parallel && is.null(n.cores)){
        n.cores <- detectCores()
        n.cores <- min(n.cores, 6)
    }


    ## GET FILE NAMES ##
    res.file.names <- paste("run", 1:n.runs, "-", res.file.name, sep="")
    tune.file.names <- paste("run", 1:n.runs, "-", tune.file.name, sep="")


    ## HANDLE SEED ##
    if(is.null(seed)){
        seed <- as.integer(runif(n.runs,min=0,max=2e9))
    } else {
        seed <- rep(seed, length=n.runs)
    }


    ## COMPUTATIONS ##
    if(parallel){
        ## create cluster ##
        clust <- makeCluster(n.cores)

        ## load outbreaker for each child ##
        clusterEvalQ(clust, library(outbreaker))

        ## transfer data onto each child ##
        listArgs <- c("dna", "dates", "idx.dna", "mut.model", "spa.model", "w.dens", "w.trunc", "f.dens", "f.trunc", "dist.mat", "init.tree", "init.kappa", "n.iter",
                      "sample.every", "tune.every", "burnin", "import.method", "find.import.n", "pi.prior1", "pi.prior2", "init.mu1", "init.mu2",
                      "init.spa1", "move.mut", "spa1.prior", "move.mut", "move.ances", "move.kappa", "move.Tinf", "move.pi", "move.spa",
                      "outlier.threshold", "max.kappa", "res.file.names", "tune.file.names", "seed")

        clusterExport(clust, listArgs, envir=environment())

        ## set calls to outbreaker on each child ##
        res <- parLapply(clust, 1:n.runs, function(i)  outbreaker(dna=dna, dates=dates, idx.dna=idx.dna,
                                                                  mut.model=mut.model, spa.model=spa.model,
                                                                  w.dens=w.dens, w.trunc=w.trunc,
                                                                  f.dens=f.dens, f.trunc=f.trunc,
                                                                  dist.mat=dist.mat, ## locations=locations,
                                                                  init.tree=init.tree, init.kappa=init.kappa,
                                                                  n.iter=n.iter, sample.every=sample.every,
                                                                  tune.every=tune.every, burnin=burnin,
                                                                  import.method=import.method,
                                                                  find.import.n=find.import.n,
                                                                  pi.prior1=pi.prior1, pi.prior2=pi.prior2,
                                                                  spa1.prior=spa1.prior,
                                                                  init.mu1=init.mu1, init.mu2=init.mu2, init.spa1=init.spa1,
                                                                  move.mut=move.mut, move.ances=move.ances, move.kappa=move.kappa,
                                                                  move.Tinf=move.Tinf, move.pi=move.pi, move.spa=move.spa,
                                                                  outlier.threshold = outlier.threshold, max.kappa=max.kappa,
                                                                  quiet=TRUE, res.file.name=res.file.names[i],
                                                                  tune.file.name=tune.file.names[i], seed=seed[i]))

        ## close parallel processes ##
        stopCluster(clust)

        ## Version with mclapply - doesn't work on windows ##
        ## res <- mclapply(1:n.runs, function(i)  outbreaker(dna=dna, dates=dates, idx.dna=idx.dna, w.dens=w.dens, w.trunc=w.trunc,
        ##                                                     init.tree=init.tree, init.kappa=init.kappa,
        ##                                                     n.iter=n.iter, sample.every=sample.every,
        ##                                                     tune.every=tune.every, burnin=burnin,
        ##                                                     find.import=find.import, find.import.n=find.import.n,
        ##                                                     pi.prior1=pi.prior1, pi.prior2=pi.prior2,
        ##                                                     init.mu1=init.mu1, init.mu2=init.mu2,
        ##                                                     move.mut=move.mut, move.ances=move.ances, move.kappa=move.kappa,
        ##                                                     move.Tinf=move.Tinf, move.pi=move.pi,
        ##                                                     quiet=TRUE, res.file.name=res.file.names[i],
        ##                                                     tune.file.name=tune.file.names[i], seed=seed[i]),
        ##                   mc.cores=n.cores, mc.silent=FALSE, mc.cleanup=TRUE, mc.preschedule=TRUE, mc.set.seed=TRUE)
    } else {
        res <- lapply(1:n.runs, function(i)  outbreaker(dna=dna, dates=dates, idx.dna=idx.dna,
                                                        mut.model=mut.model, spa.model=spa.model,
                                                        w.dens=w.dens, w.trunc=w.trunc,
                                                        f.dens=f.dens, f.trunc=f.trunc,
                                                        dist.mat=dist.mat,
                                                        init.tree=init.tree, init.kappa=init.kappa,
                                                        n.iter=n.iter, sample.every=sample.every,
                                                        tune.every=tune.every, burnin=burnin,
                                                        import.method=import.method,
                                                        find.import.n=find.import.n,
                                                        pi.prior1=pi.prior1, pi.prior2=pi.prior2,
                                                        spa1.prior=spa1.prior,
                                                        init.mu1=init.mu1, init.mu2=init.mu2, init.spa1=init.spa1,
                                                        move.mut=move.mut, move.ances=move.ances, move.kappa=move.kappa,
                                                        move.Tinf=move.Tinf, move.pi=move.pi, move.spa=move.spa,
                                                        outlier.threshold = outlier.threshold, max.kappa=max.kappa,
                                                        quiet=TRUE, res.file.name=res.file.names[i],
                                                        tune.file.name=tune.file.names[i], seed=seed[i]))
    }


    ## MERGE RESULTS ##
    res.old <- res
    res <- res[[1]]
    res$tune.end <- max(sapply(res.old, function(e) e$tune.end))
    res$chains <- Reduce(rbind, lapply(res.old, function(e) e$chains))
    res$chains$run <- factor(rep(1:n.runs, each=nrow(res.old[[1]]$chains)))
    res$n.runs <- n.runs
    res$call <- match.call()

    ## RETURN RESULTS ##
    return(res)
} # end outbreaker.parallel

