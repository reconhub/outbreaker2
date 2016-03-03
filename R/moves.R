#' Movements of augmented data and parameters for outbreaker2
#'
#' These functions are used to move in the parameter space in the MCMC.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @rdname moves
#'
#' @param data a list of data items as returned by \code{outbreaker.data}
#' @param param a list of parameters as returned by \code{outbreaker.mcmc.init}
#' @param rand  a list of items as returned by \code{outbreaker.rand.vec}
#'
#' @return a potentially modified list of parameters as returned by \code{outbreaker.mcmc.init}
#'
#' @importFrom stats rnorm
#'
move.mu <- function(data, param, config, rand){
    ## get new proposed values
    new.param <- param
    new.param$current.mu <- new.param$current.mu + rand$mu.rnorm1()

    ## escape if new.mu<0 or >1
    if(new.param$current.mu<0 || new.param$current.mu>1) return(param)

    ## compute log ratio  (assumes symmetric proposal)
    logratio <- post.genetic(data=data, param=new.param) -
        post.genetic(data=data, param=param)

    ## accept/reject
    if(logratio >= rand$log.runif1()) return(new.param)
    return(param)
} # end move.mu





#' @rdname moves
#' @export
#'
move.t.inf <- function(data, param, config, rand){ # assumes symmetric proposal

    ## propose new t.inf
    new.param <- param
    new.param$current.t.inf <- new.param$current.t.inf +
        sample(-1:1, size=length(param$current.t.inf), replace=TRUE, prob=c(.1,8,.1))

    ## compute log ratio
    logratio <- ll.timing(data=data, param=new.param) - ll.timing(data=data, param=param)

    ## accept/reject
    if(logratio >= rand$log.runif1()) return(new.param)

    return(param)
} # end move.t.inf




#' @rdname moves
#' @export
#' @param config a list of settings as returned by \code{outbreaker.config}
#'
move.ances <- function(data, param, config, rand){
    ## create new parameters
    new.param <- param

    ## find out which ancestries to move
    ances.can.move <- !is.na(param$current.ances) & param$current.t.inf>min(param$current.t.inf)
    if(!any(ances.can.move)){
        warning("trying to move ancestries but none can move")
        return(param$current.ances)
    }
    n.to.move <- max(round(config$prop.ances.move * sum(ances.can.move)),1)
    to.move <- sample(which(ances.can.move), n.to.move, replace=FALSE)

    ## initialize new ances
    new.param$current.ances <- param$current.ances

    ## move all ancestries that should be moved
    for(i in to.move){
        ## propose new ancestor
        new.param$current.ances[i] <- find.possible.ances(param$current.t.inf, i)

        ## compute log ratio
        logratio <-  ll.all(data=data, param=new.param) - ll.all(data=data, param=param)

        ## accept/reject
        if(logratio >= rand$log.runif1()){
            param$current.ances[i] <- new.param$current.ances[i]
        } else {
            new.param$current.ances[i] <- param$current.ances[i]
        }
    } # end for loop

    return(param)
} # end move.ances






#' @rdname moves
#' @export
#'
move.swap.ances <- function(data, param, config, rand){
    ## create new parameters
    new.param <- param

    ## find out which ancestries to move
    ances.can.move <- !is.na(param$current.ances) & param$current.t.inf>min(param$current.t.inf)
    if(!any(ances.can.move)){
        warning("trying to move ancestries but none can move")
        return(param$current.ances)
    }
    n.to.move <- max(round(config$prop.ances.move * sum(ances.can.move)),1)
    to.move <- sample(which(ances.can.move), n.to.move, replace=FALSE)

    ## initialize new ances and t.inf
    new.param$current.ances <- param$current.ances
    new.param$current.t.inf <- param$current.t.inf

    ## move all ancestries that should be moved
    for(i in to.move){
        ## swap ancestries
        new.param <- swap.ances(new.param, i)

        ## compute log ratio
        logratio <- ll.all(data=data, param=new.param) - ll.all(data=data, param=param)

        ## accept/reject
        if(logratio >= rand$log.runif1()){
            param$current.ances[i] <- new.param$current.ances[i]
            param$current.t.inf[i] <- new.param$current.t.inf[i]
        } else {
            new.param$current.ances[i] <- param$current.ances[i]
            new.param$current.t.inf[i] <- param$current.t.inf[i]
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
swap.ances <- function(param, i){
    ## stop if 'i' out of range
    if(i>length(param$current.ances)) stop("trying to swap ancestry of case ", i, " while there are only ", length(param$current.ances), " cases")

    ## find ancestor of 'i'
    x <- param$current.ances[i]

    ## stop if case 'i' is imported
    if(is.na(x)){
        warning("trying to swap the ancestry of the imported case ", i)
        return(param)
    }

    ## find indices to swap
    to.be.x <- which(param$current.ances==i)
    to.be.i <- which(param$current.ances==x)

    ## swap 'i' and 'x' in ancestries
    param$current.ances[to.be.x] <- x
    param$current.ances[to.be.i] <- i

    ## the ancestor of 'i' is now has the ancestor of 'x'
    param$current.ances[i] <- param$current.ances[x]

    ## 'i' is now the ancestor of 'x'
    param$current.ances[x] <- i

    ## swap t.inf
    param$current.t.inf[c(x,i)] <- param$current.t.inf[c(i,x)]

    ## return
    return(param)
} # end swap.ancestries
