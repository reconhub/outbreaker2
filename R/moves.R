
## In all of these functions:
## --------------------------
##
## We return a function in which all elements required for moving parameters are enclosed; this
## includes:

## - config: a list containing info about general settings of the method

## - densities: a named list with 3 components (loglike, priors, posteriors), each being a list of
## functions with a single argument 'param' (the rest is enclosed in the functions

## - rand: a list containing pre-generated random numbers

## See 'likelihood.R' for more detail about the motivation (but basically, it's faster). These
## functions are later called by 'create.moves()' which will make a list of functions with enclosed
## items.




## Movement of the mutation rate 'mu' is done using a dumb normal proposal. This is satisfying for
## now - we only reject a few non-sensical values outside the range [0;1]. The SD of the proposal
## (implicitely contained in rand$mu.rnorm1, but really provided through 'config', seems fine as the
## range of real values will never change much. Probably not much point in using auto-tuning here.

make.move.mu <- function(config, densities) {
    data <- environment(densities$loglike$genetic)$data
    ## .move.mu(config, densities) # uncomment for pure R version
    function(param) {
        cpp.move.mu(data, param, config)
        return(param)
    }
}




## Movement of infection dates are +/- 1 from current states. These movements are currently
## vectorised, i.e. a bunch of dates are proposed all together; this may not be sustainable for
## larger datasets. The non-vectorised option will be slower and speed-up with C/C++ will be more
## substantial then.

make.move.t.inf <- function(config, densities) {
    data <- environment(densities$loglike$timing)$data
    function(param) {
        cpp.move.t.inf(data, param)
        return(param)
    }
}





## Movement of ancestries ('alpha') is not vectorised, movements are made one case at a time. This
## procedure is simply about picking an infector at random amongst cases preceeding the case
## considered. This movement is not symmetric, as the number of choices may change. The original
## version in 'outbreaker' used to move simultaneously 'alpha', 'kappa' and 't.inf', but current
## implementation is simpler and seems to mix at least as well. Proper movement of 'alpha' needs
## this procedure as well as a swapping procedure (swaps are not possible through move.alpha only).

make.move.alpha <- function(config, densities) {
    function(param) {
        ## create new parameters
        new.param <- param

        ## find out which ancestries to move
        alpha.can.move <- !is.na(param$current.alpha) & param$current.t.inf>min(param$current.t.inf)
        if (!any(alpha.can.move)) {
            warning("trying to move ancestries but none can move")
            return(param$current.alpha)
        }
        n.to.move <- max(round(config$prop.alpha.move * sum(alpha.can.move)),1)
        to.move <- sample(which(alpha.can.move), n.to.move, replace=FALSE)

        ## initialize new alpha
        new.param$current.alpha <- param$current.alpha

        ## move all ancestries that should be moved
        for (i in to.move) {
            ## propose new ancestor
            new.param$current.alpha[i] <- .choose.possible.alpha(param$current.t.inf, i)

            ## compute log ratio
            logratio <-  densities$loglike$all(new.param) - densities$loglike$all(param)

            ## compute correction factor
            logratio <- logratio + log(sum(.are.possible.alpha(new.param$current.t.inf, i))) -
                log(sum(.are.possible.alpha(param$current.t.inf, i)))

            ## accept/reject
            if (logratio >= log(stats::runif(1))) {
                param$current.alpha[i] <- new.param$current.alpha[i]
            } else {
                new.param$current.alpha[i] <- param$current.alpha[i]
            }
        } # end for loop

        return(param)
    }
}





## This is the complementary procedure to the above one (move.alpha). This type of move swaps a case
## 'a' with its ancestor, e.g.

## x -> a -> b  becomes a -> x -> b

## Obviously cases are moved one at a time. We need to used local likelihood changes for this move
## to scale well with outbreak size. The complicated bit is that the move impacts all descendents
## from 'a' as well as 'x'.

make.move.swap.cases <- function(config, densities) {
    function(param) {
        ## find ancestries which can move
        to.move <- select.alpha.to.move(param, config)

        ## leave if nothing moves
        if (length(to.move)<1) {
            return(param)
        }

        ## move all ancestries that should be moved
        for (i in to.move) {
            ## swap ancestries
            new.param <- swap.cases(param, config, i)

            ## compute log ratio using local changes only; these include:

            ## descendents of to.move
            ## descendents of alpha[to.move]
            ## alpha[to.move]

            affected.cases <- c(find.descendents(param, i=i),
                                find.descendents(param, i=param$current.alpha[i]),
                                param$current.alpha[i])
            logratio <- densities$loglike$all(new.param, i=affected.cases) - densities$loglike$all(param, i=affected.cases)

            ## accept/reject
            if (logratio >= log(stats::runif(1))) {
                param <- new.param
            }
        } # end for loop

        return(param)
    }
}




## This movement of the reporting probability 'pi' is treated in a similar way to the mutation rate
## 'mu'; the only difference lies in the default SD for the proposal Normal distribution, larger for
## 'pi' (jumps are expected to be a bit larger as values are typically larger for 'pi'). Again, not
## many dumb values proposed here, so it may not be worth it to use a non-symmetric proposal
## (e.g. log-normal).

make.move.pi <- function(config, densities) {
    function(param) {
        ## get new proposed values
        new.param <- param
        ## new.param$current.pi <- new.param$current.pi + rand$pi.rnorm1()
        new.param$current.pi <- stats::rnorm(1, mean=new.param$current.pi, sd=config$sd.pi)

        ## escape if new.pi<0 or >1
        if (new.param$current.pi<0 || new.param$current.pi>1) {
            return(param)
        }

        ## compute log ratio  (assumes symmetric proposal)
        logratio <- densities$posteriors$reporting(new.param) -
            densities$posteriors$reporting(param)

        ## accept/reject
        if (logratio >= log(stats::runif(1))) {
            return(new.param)
        }
        return(param)
    }
}





## Movement of the number of generations on transmission chains ('kappa') is done for one ancestry
## at a time. As for infection times ('t.inf') we use a dumb, symmetric +/- 1 proposal. But because
## values are typically in a short range (e.g. [1-3]) we probably propose more dumb values here. We
## may eventually want to bounce back or use and correct for assymetric proposals.

make.move.kappa <- function(config, densities) {
    function(param) {

        ## determine which cases to move
        kappa.can.move <- !is.na(param$current.kappa)
        n.to.move <- max(round(.2 * sum(kappa.can.move), 1))
        to.move <- sample(which(kappa.can.move), n.to.move, replace=FALSE)

        ## initialize new kappa
        new.param <- param

        ## move all ancestries that should be moved
        for (i in to.move) {
            ## propose new kappa
            new.param$current.kappa[i] <- new.param$current.kappa[i] + sample(c(-1,1), size=1)

            ## reject move automatically if new kappa < 1 or greater than allowed max
            if (new.param$current.kappa[i] < 1 ||
               new.param$current.kappa[i] > config$max.kappa) {
                new.param$current.kappa[i] <- param$current.kappa[i]
            } else {
                ## compute log ratio
                logratio <- densities$loglike$timing.infections(new.param, i) +
                    densities$loglike$genetic(new.param, i) +
                        densities$loglike$reporting(new.param, i) -
                            densities$loglike$timing.infections(param, i) -
                                densities$loglike$genetic(param, i) -
                                    densities$loglike$reporting(param, i)

                ## accept/reject
                if (logratio >= log(stats::runif(1))) {
                    param$current.kappa[i] <- new.param$current.kappa[i]
                } else {
                    new.param$current.kappa[i] <- param$current.kappa[i]
                }
            }
        } # end for loop


        ## output is a list of potentially modified parameters
        return(param)
    }
}
