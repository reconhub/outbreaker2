## #' Random number pre-generation for outbreaker2
## #'
## #' This function generates pre-drawn vectors of random numbers to be used in outbreaker2.
## #'
## #' @author Thibaut Jombart \email{t_jombart@@imperial_ac_uk}
## #'
## #' @rdname outbreaker_rand_vec
## #'
## #' @param config a list of settings as returned by \code{outbreaker_config}
## #'
## #' @export
## #'
## outbreaker_rand_vec <- function(config)
## {
##     ## CREATE OUTPUT ##
##     out <- list()

##     ## CREATE FAST RANDOM NUMBER GENERATORS ##
##     out$log_runif1 <- make_fast_rand1(f = runif, batch_size = config$batch_size, log = TRUE)
##     out$mu_rnorm1 <- make_fast_rand1(mean = 0, sd = config$sd_mu, f = rnorm, batch_size = config$batch_size, log = FALSE)
##     out$pi_rnorm1 <- make_fast_rand1(mean = 0, sd = config$sd_pi, f = rnorm, batch_size = config$batch_size, log = FALSE)

##     return(out)
## }
