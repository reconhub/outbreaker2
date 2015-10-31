#' Movements of augmented data and parameters for outbreaker2
#'
#' These functions are used to move in the parameter space in the MCMC.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @rdname moves
#'
#' @param data a list of data items as returned by \code{outbreaker.data}
#' @param chain a list of output items as returned by \code{outbreaker.mcmc.init}
#' @param r.new a deviation from the previous mu
#' @param r.acc a logged random variate from Unif(0,1)
#'
#' @importFrom stats rnorm
#'
move.mu <- function(data, chain, r.new, r.acc){
    ## get new proposed values
    new.mu <- chain$current.mu + r.new

    ## escape if new.mu<0 or >1
    if(new.mu<0 || new.mu>1) return(chain$current.mu)

    ## compute log ratio  (assumes symmetric proposal)
    logratio <- post.genetic(D=data$D, gen.length=data$L, ances=chain$current.ances, mu=new.mu) -
        post.genetic(D=data$D, gen.length=data$L, ances=chain$current.ances, mu=chain$current.mu)

    ## accept/reject
    if(logratio >= r.acc) return(new.mu)
    return(chain$current.mu)
} # end move.mu





#' @rdname moves
#' @export
#'
move.t.inf <- function(data, chain, r.acc){ # assumes symmetric proposal
    ## propose new t.inf
    new.t.inf <- chain$current.t.inf + sample(-1:1, size=length(chain$current.t.inf), replace=TRUE, prob=c(.1,8,.1))

    ## compute log ratio
    logratio <- ll.timing(log.w=data$log.w, log.f=data$log.f, sampling.times=data$sampling.times,
                          ances=chain$ances, t.inf=new.t.inf) -
        ll.timing(log.w=data$log.w, log.f=data$log.f, sampling.times=data$sampling.times,
                          ances=chain$ances, t.inf=chain$current.t.inf)

    ## accept/reject
    if(logratio >= r.acc) return(new.t.inf)

    return(chain$current.t.inf)
} # end move.t.inf





#' @rdname moves
#' @export
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





#' @rdname moves
#' @export
#'
move.ances <- function(t.inf, sampling.times, D, gen.length, log.w, log.f, ances, mu, lunif){
    ## propose new ances
    new.ances <- rances(t.inf)

    ## compute log ratio
    logratio <- ll.all(t.inf=t.inf, sampling.times=sampling.times, D=D, gen.length=gen.length,
                       log.w=log.w, log.f=log.f, ances=new.ances, mu=mu) -
                           ll.all(t.inf=t.inf, sampling.times=sampling.times, D=D, gen.length=gen.length,
                                  log.w=log.w, log.f=log.f, ances=ances, mu=mu)

    ## accept/reject
    if(logratio >= lunif) return(new.ances)

    return(ances)
} # end move.ances
