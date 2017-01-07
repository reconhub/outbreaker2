
## In all of these functions:
## --------------------------
##
## We return a function in which all elements required for moving parameters are
## enclosed; this includes:

## - config: a list containing info about general settings of the method

## - densities: a named list with 3 components (loglike, priors, posteriors),
## each being a list of functions with a single argument 'param' (the rest is
## enclosed in the functions

## - rand: a list containing pre-generated random numbers

## See 'likelihood.R' for more detail about the motivation (but basically, it's
## faster). These functions are later called by 'create.moves()' which will make
## a list of functions with enclosed items.




## Movement of the mutation rate 'mu' is done using a dumb normal proposal. This
## is satisfying for now - we only reject a few non-sensical values outside the
## range [0;1]. The SD of the proposal (implicitely contained in rand$mu.rnorm1,
## but really provided through 'config', seems fine as the range of real values
## will never change much. Probably not much point in using auto-tuning here.

make.move.mu <- function(config, densities) {
    data <- environment(densities$loglike$genetic)$data
    ## .move.mu(config, densities) # uncomment for pure R version
    function(param) {
        cpp_move_mu(data, param, config)
    }
}




## Movement of infection dates are +/- 1 from current states. These movements
## are currently vectorised, i.e. a bunch of dates are proposed all together;
## this may not be sustainable for larger datasets. The non-vectorised option
## will be slower and speed-up with C/C++ will be more substantial then.

make.move.t.inf <- function(config, densities) {
    data <- environment(densities$loglike$timing)$data
    ## .move.t.inf(config, densities)
    function(param) {
        cpp_move_t_inf(data, param)
    }
}





## Movement of ancestries ('alpha') is not vectorised, movements are made one
## case at a time. This procedure is simply about picking an infector at random
## amongst cases preceeding the case considered. This movement is not symmetric,
## as the number of choices may change. The original version in 'outbreaker'
## used to move simultaneously 'alpha', 'kappa' and 't.inf', but current
## implementation is simpler and seems to mix at least as well. Proper movement
## of 'alpha' needs this procedure as well as a swapping procedure (swaps are
## not possible through move.alpha only).

make.move.alpha <- function(config, densities) {
    data <- environment(densities$loglike$timing)$data
    ## .move.alpha(config, densities)
    function(param) {
        cpp_move_alpha(data, param)
    }
}






## This is the complementary procedure to the above one (move.alpha). This type
## of move swaps a case 'a' with its ancestor, e.g.

## x -> a -> b  becomes a -> x -> b

## Obviously cases are moved one at a time. We need to used local likelihood
## changes for this move to scale well with outbreak size. The complicated bit
## is that the move impacts all descendents from 'a' as well as 'x'.

make.move.swap.cases <- function(config, densities) {
    data <- environment(densities$loglike$timing)$data
    ##.move.swap.cases(config, densities) # uncomment for pure R version
    function(param) {
        cpp_move_swap_cases(data, param)
    }
}





## This movement of the reporting probability 'pi' is treated in a similar way
## to the mutation rate 'mu'; the only difference lies in the default SD for the
## proposal Normal distribution, larger for 'pi' (jumps are expected to be a bit
## larger as values are typically larger for 'pi'). Again, not many dumb values
## proposed here, so it may not be worth it to use a non-symmetric proposal
## (e.g. log-normal).

make.move.pi <- function(config, densities) {
    data <- environment(densities$loglike$genetic)$data
    ## .move.pi(config, densities) # uncomment for pure R version
    function(param) {
        cpp_move_pi(data, param, config)
    }

}





## Movement of the number of generations on transmission chains ('kappa') is
## done for one ancestry at a time. As for infection times ('t.inf') we use a
## dumb, symmetric +/- 1 proposal. But because values are typically in a short
## range (e.g. [1-3]) we probably propose more dumb values here. We may
## eventually want to bounce back or use and correct for assymetric proposals.

make.move.kappa <- function(config, densities) {
    data <- environment(densities$loglike$genetic)$data
    .move.kappa(config, densities) # uncomment for pure R version
    ## function(param) {
    ##     cpp_move_kappa(data, param, config)
    ## }

}
