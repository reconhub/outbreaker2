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
move.mu <- function(param, config, densities, rand){
    ## get new proposed values
    new.param <- param
    new.param$current.mu <- new.param$current.mu + rand$mu.rnorm1()

    ## escape if new.mu<0 or >1
    if(new.param$current.mu<0 || new.param$current.mu>1) return(param)

    ## compute log ratio  (assumes symmetric proposal)
    logratio <- densities$posteriors$post.genetic(new.param) -
        densities$posteriors$post.genetic(param)

    ## accept/reject
    if(logratio >= rand$log.runif1()) return(new.param)
    return(param)
} # end move.mu





#' @rdname moves
#' @export
#'
move.t.inf <- function(param, config, densities, rand){ # assumes symmetric proposal

    ## propose new t.inf
    new.param <- param
    new.param$current.t.inf <- new.param$current.t.inf +
        sample(-1:1, size=data$N, replace=TRUE, prob=c(.1,8,.1))

    ## compute log ratio
    logratio <- densities$loglike$ll.timing(new.param) - densities$loglike$ll.timing(param)

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
move.alpha <- function(param, config, densities, rand){
    ## create new parameters
    new.param <- param

    ## find out which ancestries to move
    alpha.can.move <- !is.na(param$current.alpha) & param$current.t.inf>min(param$current.t.inf)
    if(!any(alpha.can.move)){
        warning("trying to move ancestries but none can move")
        return(param$current.alpha)
    }
    n.to.move <- max(round(config$prop.alpha.move * sum(alpha.can.move)),1)
    to.move <- sample(which(alpha.can.move), n.to.move, replace=FALSE)

    ## initialize new alpha
    new.param$current.alpha <- param$current.alpha

    ## move all ancestries that should be moved
    for(i in to.move){
        ## propose new ancestor
        new.param$current.alpha[i] <- choose.possible.alpha(param$current.t.inf, i)

        ## compute log ratio
        logratio <-  densities$loglike$ll.all(new.param) - densities$loglike$ll.all(param)

        ## compute correction factor
        logratio <- logratio + log(sum(are.possible.alpha(new.param$current.t.inf, i))) -
            log(sum(are.possible.alpha(param$current.t.inf, i)))

        ## accept/reject
        if(logratio >= rand$log.runif1()){
            param$current.alpha[i] <- new.param$current.alpha[i]
        } else {
            new.param$current.alpha[i] <- param$current.alpha[i]
        }
    } # end for loop

    return(param)
} # end move.alpha






#' @rdname moves
#' @export
#'
move.swap.cases <- function(param, config, densities, rand){
    ## find ancestries which can move
    to.move <- select.alpha.to.move(param, config)

    ## leave if nothing moves
    if(length(to.move)<1) return(param)

    ## move all ancestries that should be moved
    for(i in to.move){
        ## swap ancestries
        new.param <- swap.cases(param, config, i)

        ## compute log ratio
        ## only use local changes:
        ## descendents of to.move
        ## descendents of alpha[to.move]
        ## alpha[to.move]
        affected.cases <- c(find.descendents(param, i=i),
                            find.descendents(param, i=param$current.alpha[i]),
                            param$current.alpha[i])
        logratio <- densities$loglike$ll.all(new.param, i=affected.cases) - densities$loglike$ll.all(param, i=affected.cases)

        ## accept/reject
        if(logratio >= rand$log.runif1()){
            param <- new.param
        }
    } # end for loop

    return(param)
} # end move.swap.cases





#' @rdname moves
#' @export
#'
move.pi <- function(param, config, densities, rand){
    ## get new proposed values
    new.param <- param
    new.param$current.pi <- new.param$current.pi + rand$pi.rnorm1()

    ## escape if new.pi<0 or >1
    if(new.param$current.pi<0 || new.param$current.pi>1) return(param)

    ## compute log ratio  (assumes symmetric proposal)
    logratio <- densities$posteriors$post.reporting(new.param) -
        densities$posteriors$post.reporting(param)

    ## accept/reject
    if(logratio >= rand$log.runif1()) return(new.param)
    return(param)
} # end move.pi





#' @rdname moves
#' @export
#'
move.kappa <- function(param, config, densities, rand){

    ## determine which cases to move
    kappa.can.move <- !is.na(param$current.kappa)
    n.to.move <- max(round(.2 * sum(kappa.can.move), 1))
    to.move <- sample(which(kappa.can.move), n.to.move, replace=FALSE)

    ## initialize new kappa
    new.param <- param

    ## move all ancestries that should be moved
    for(i in to.move){
        ## propose new kappa
        new.param$current.kappa[i] <- new.param$current.kappa[i] + sample(c(-1,1), size=1)

        ## reject move automatically if new kappa < 1 or greater than allowed max
        if(new.param$current.kappa[i] < 1 ||
           new.param$current.kappa[i] > config$max.kappa){
            new.param$current.kappa[i] <- param$current.kappa[i]
        } else {
            ## compute log ratio
            logratio <- densities$loglike$ll.timing.infections(new.param, i) +
                densities$loglike$ll.genetic(new.param, i) +
                densities$loglike$ll.reporting(new.param, i) -
                densities$loglike$ll.timing.infections(param, i) -
                densities$loglike$ll.genetic(param, i) -
                densities$loglike$ll.reporting(param, i)

            ## accept/reject
            if(logratio >= rand$log.runif1()){
                param$current.kappa[i] <- new.param$current.kappa[i]
            } else {
                new.param$current.kappa[i] <- param$current.kappa[i]
            }
        }
    } # end for loop


    ## return potentially modified parameters
    return(param)
} # end move.kappa
