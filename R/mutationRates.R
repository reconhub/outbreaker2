############
## get.mu ##
###########
get.mu <- function(x, burnin=2e4, genome.size=NULL){

    ## GET DATA ##
    ## remove burnin
    if(!any(x$chains$step>burnin)) stop("burnin too high - no chain left")
    dat <- x$chains[x$chains$step>burnin,,drop=FALSE]
    ances <- dat[,grep("alpha", names(x$chains)),drop=FALSE] # table of ancestries
    Tinf <-  dat[,grep("Tinf", names(x$chains)),drop=FALSE]
    D <- as.matrix(x$D)
    n <- ncol(ances)


    ## AUXILIARY FUNCTIONS, GET RESULTS FOR ONE SAMPLE ##
    f1 <- function(vecAnces, vecTinf){
        ## get the number of mutations between cases ##
        vecAnces <- as.integer(vecAnces)
        vecAnces[vecAnces<1] <- NA
        nmut <- sapply(1:n, function(i) D[i,vecAnces[i]])

        ## get the time between cases ##
        vecTinf <- as.integer(vecTinf)
        deltaT <- vecTinf-vecTinf[vecAnces]

        ## get mutation rate ##
        out <- mean(nmut/deltaT, na.rm=TRUE)
        return(out)
    }

    ## GET MUTATION RATES FROM POSTERIOR SAMPLES ##
    out <- unlist(lapply(1:nrow(ances), function(i) f1(ances[i,],Tinf[i,])))

    ## rescale mutation rate if necessary ##
    if(!is.null(genome.size)) out <- out/genome.size
    return(out)
} # end get.mu
