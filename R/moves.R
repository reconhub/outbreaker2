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
        sample(-1:1, size=data$N, replace=TRUE, prob=c(.1,8,.1))

    ## compute log ratio
    logratio <- ll.timing(data=data, param=new.param) - ll.timing(data=data, param=param)

    ## accept/reject
    if(logratio >= rand$log.runif1()){
        return(new.param)
    } else {
        return(param)
    }
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
        new.param$current.ances[i] <- choose.possible.ances(param$current.t.inf, i)

        ## compute log ratio
        logratio <-  ll.all(data=data, param=new.param) - ll.all(data=data, param=param)

        ## compute correction factor
        logratio <- logratio + log(sum(are.possible.ances(new.param$current.t.inf, i))) -
            log(sum(are.possible.ances(param$current.t.inf, i)))

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
    ## find ancestries which can move
    to.move <- select.ances.to.move(param, config)

    ## leave if nothing moves
    if(length(to.move)<1) return(param)

    ## move all ancestries that should be moved
    for(i in to.move){
        ## swap ancestries
        new.param <- swap.ances(param, config, i)

        ## compute log ratio
        ## only use local changes:
        ## descendents of to.move
        ## descendents of ances[to.move]
        ## ances[to.move]
        affected.cases <- c(find.descendents(data=data, param=param, i=i),
                            find.descendents(data=data, param=param, i=param$current.ances[i]),
                            param$current.ances[i])
        logratio <- ll.all(data=data, param=new.param, i=affected.cases) - ll.all(data=data, param=param, i=affected.cases)

        ## accept/reject
        if(logratio >= rand$log.runif1()){
            param <- new.param
        }
    } # end for loop

    return(param)
} # end move.swap.ances





#' @rdname moves
#' @export
#' @param t.inf a vector of infection dates
#'
rances <- function(t.inf){
    ## choose.possible.ancestors
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
move.pi <- function(data, param, config, rand){
    ## get new proposed values
    new.param <- param
    new.param$current.pi <- new.param$current.pi + rand$pi.rnorm1()

    ## escape if new.pi<0 or >1
    if(new.param$current.pi<0 || new.param$current.pi>1) return(param)

    ## compute log ratio  (assumes symmetric proposal)
    logratio <- post.genetic(data=data, param=new.param) -
        post.genetic(data=data, param=param)

    ## accept/reject
    if(logratio >= rand$log.runif1()) return(new.param)
    return(param)
} # end move.pi

