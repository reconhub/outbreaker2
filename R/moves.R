#' Movements of augmented data and parameters for outbreaker2
#'
#' These functions are used to move in the parameter space in the MCMC.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @rdname moves
#'
#' @param data a list of data items as returned by \code{outbreaker.data}
#' @param param a list of output items as returned by \code{outbreaker.mcmc.init}
#' @param rand  a list of items as returned by \code{outbreaker.rand.vec}
#' @importFrom stats rnorm
#'
move.mu <- function(data, param, rand){
    ## get new proposed values
    new.mu <- param$current.mu + rand$mu.rnorm1()

    ## escape if new.mu<0 or >1
    if(new.mu<0 || new.mu>1) return(param$current.mu)

    ## compute log ratio  (assumes symmetric proposal)
    logratio <- post.genetic(D=data$D, gen.length=data$L, ances=param$current.ances, mu=new.mu) -
        post.genetic(D=data$D, gen.length=data$L, ances=param$current.ances, mu=param$current.mu)

    ## accept/reject
    if(logratio >= rand$log.runif1()) return(new.mu)
    return(param$current.mu)
} # end move.mu





#' @rdname moves
#' @export
#'
move.t.inf <- function(data, param, rand){ # assumes symmetric proposal

    ## propose new t.inf
    new.t.inf <- param$current.t.inf + sample(-1:1, size=length(param$current.t.inf), replace=TRUE, prob=c(.1,8,.1))

    ## compute log ratio
    logratio <- ll.timing(log.w=data$log.w, log.f=data$log.f, sampling.times=data$sampling.times,
                          ances=param$current.ances, t.inf=new.t.inf) -
                              ll.timing(log.w=data$log.w, log.f=data$log.f, sampling.times=data$sampling.times,
                                        ances=param$current.ances, t.inf=param$current.t.inf)

    ## accept/reject
    if(logratio >= rand$log.runif1()) return(new.t.inf)

    return(param$current.t.inf)
} # end move.t.inf




#' @rdname moves
#' @export
#' @param config a list of settings as returned by \code{outbreaker.config}
move.ances <- function(data, param, config, rand){
    ## find out which ancestries to move
    ances.can.move <- !is.na(param$current.ances) & param$current.t.inf>min(param$current.t.inf)
    if(!any(ances.can.move)){
        warning("trying to move ancestries but none can move")
        return(param$current.ances)
    }
    n.to.move <- max(round(config$prop.ances.move * sum(ances.can.move)),1)
    to.move <- sample(which(ances.can.move), n.to.move, replace=FALSE)

    ## initialize new ances
    new.ances <- param$current.ances

    ## move all ancestries that should be moved
    for(i in to.move){
        ## propose new ancestor
        new.ances[i] <- find.possible.ances(param$current.t.inf, i)

        ## compute log ratio
        logratio <- ll.timing(log.w=data$log.w, log.f=data$log.f, sampling.times=data$sampling.times,
                              t.inf=param$current.t.inf, ances=new.ances) +
                                  ll.genetic(D=data$D, gen.length=data$L, mu=param$current.mu, ances=new.ances) -
                                      ll.timing(log.w=data$log.w, log.f=data$log.f, sampling.times=data$sampling.times,
                                                t.inf=param$current.t.inf, ances=param$current.ances) -
                                                    ll.genetic(D=data$D, gen.length=data$L, mu=param$current.mu, ances=param$current.ances)

        ## accept/reject
        if(logratio >= rand$log.runif1()){
            param$current.ances[i] <- new.ances[i]
        } else {
            new.ances[i] <- param$current.ances[i]
        }
    } # end for loop

    return(param$current.ances)
} # end move.ances






#' @rdname moves
#' @export
#'
move.swap.ances <- function(data, param, config, rand){
     ## find out which ancestries to move
    ances.can.move <- !is.na(param$current.ances) & param$current.t.inf>min(param$current.t.inf)
    if(!any(ances.can.move)){
        warning("trying to move ancestries but none can move")
        return(param$current.ances)
    }
    n.to.move <- max(round(config$prop.ances.move * sum(ances.can.move)),1)
    to.move <- sample(which(ances.can.move), n.to.move, replace=FALSE)

    ## initialize new ances and t.inf
    new.ances <- param$current.ances
    new.t.inf <- param$current.t.inf

    ## move all ancestries that should be moved
    for(i in to.move){
        ## swap ancestries
        temp <- swap.ances(new.ances, new.t.inf, i)
        new.ances <- temp$ances
        new.t.inf <- temp$t.inf

        ## compute log ratio
        logratio <- ll.timing(log.w=data$log.w, log.f=data$log.f, sampling.times=data$sampling.times,
                              t.inf=new.t.inf, ances=new.ances) +
                                  ll.genetic(D=data$D, gen.length=data$L, mu=param$current.mu, ances=new.ances) -
                                      ll.timing(log.w=data$log.w, log.f=data$log.f, sampling.times=data$sampling.times,
                                                t.inf=param$current.t.inf, ances=param$current.ances) -
                                                    ll.genetic(D=data$D, gen.length=data$L, mu=param$current.mu, ances=param$current.ances)

        ## accept/reject
        if(logratio >= rand$log.runif1()){
            param$current.ances[i] <- new.ances[i]
            param$current.t.inf[i] <- new.t.inf[i]
        } else {
            new.ances[i] <- param$current.ances[i]
            new.t.inf[i] <- param$current.t.inf[i]
        }
    } # end for loop

    return(param)
} # end move.swap.ances





#' @rdname moves
#' @export
#' @param t.inf a vector of infection dates
#'
rances <- function(t.inf){
    ## find possible ancestors
    canBeAnces <- outer(t.inf,t.inf,FUN="<") # strict < is needed as we impose w(0)=0
    diag(canBeAnces) <- FALSE

    ## pick possible ancestors at random
    ances <- apply(canBeAnces, 2, function(e) ifelse(length(which(e))>0, sample(which(e),1), NA) )

    ## return
    return(ances)
} # end rances





## non-exported function
## finds possible ancestor for a case 'i'
## (any case before i)
find.possible.ances <- function(t.inf, i){
    if(length(i)>1) stop("i has a length > 1")
    if(any(t.inf[i]==min(t.inf))) return(NA)
    return(sample(which(t.inf < t.inf[i[1]]), 1))
} # end find.possible.ances





## non-exported function
## swaps ancestries in the tree
## x-> i becomes i->x
## plus all subsequent changes
swap.ances <- function(ances, t.inf, i){
    ## stop if 'i' out of range
    if(i>length(ances)) stop("trying to swap ancestry of case ", i, " while there are only ", length(ances), " cases")

    ## find ancestor of 'i'
    x <- ances[i]

    ## stop if case 'i' is imported
    if(is.na(x)){
        warning("trying to swap the ancestry of the imported case ", i)
        return(list(ances=ances, t.inf=t.inf))
    }

    ## find indices to swap
    to.be.x <- which(ances==i)
    to.be.i <- which(ances==x)

    ## swap 'i' and 'x' in ancestries
    ances[to.be.x] <- x
    ances[to.be.i] <- i

    ## the ancestor of 'i' is now has the ancestor of 'x'
    ances[i] <- ances[x]

    ## 'i' is now the ancestor of 'x'
    ances[x] <- i

    ## swap t.inf
    t.inf[c(x,i)] <- t.inf[c(i,x)]

    ## return
    return(list(ances=ances, t.inf=t.inf))
} # end swap.ancestries
