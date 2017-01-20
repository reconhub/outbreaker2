context("Test prior functions")


test_that("Priors in cpp and R references give identical results.", {
    ## skip on CRAN
    skip_on_cran()
    

    ## generate data
    param <- list(mu = 0.000123, pi = 0.789)
    config <- create_config()

    mu_r <- .prior_mu(param, config$prior_mu)
    mu_cpp <- cpp_prior_mu(param, config)

    pi_r <- .prior_pi(param, config$prior_pi[1], config$prior_pi[2])
    pi_cpp <- cpp_prior_pi(param, config)

    ## checks
    expect_equal(mu_r, mu_cpp)
    expect_equal(pi_r, pi_cpp)    
})





test_that("Sum of priors is consistent.", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    param <- list(mu = 0.000123, pi = 0.789)
    config <- create_config()

    p_mu <- cpp_prior_mu(param, config)
    p_pi <- cpp_prior_pi(param, config)
    p_all<- cpp_prior_all(param, config)
    
 
    ## checks
    expect_equal(p_mu + p_pi, p_all)

})





test_that("Prior customisation.", {

    ## skip on CRAN
    skip_on_cran()

    
    ## check default config
    config <- create_config()
    expect_equal(custom_priors(),
                     custom_priors(config))

    ## check errors
    msg <- "The following priors are not functions: mu"
    expect_error(custom_priors(mu = "Chtulhu"), msg)

    msg <- "The following priors dont' have a single argument: mu"
    expect_error(custom_priors(mu = plot), msg)


    ## custom prior parameters
    mu <- 0.00123123
    pi  <-  0.8765
    config <- create_config(prior_mu = 0.01,
                            prior_pi = c(2, 1))
    param <- list(mu = mu, pi = pi)
    p_mu <- cpp_prior_mu(param, config)
    p_pi <- cpp_prior_pi(param, config)
    expect_equal(p_mu,
                 dexp(mu, 0.01, log = TRUE))
    expect_equal(p_pi,
                 dbeta(pi, 2, 1, log = TRUE))

    
    ## custom functions
    f_mu <- function(x) { dexp(x$mu, rate = 100, log = TRUE) }
    p_mu <- cpp_prior_mu(param, config, f_mu)
    p_mu_ref <- f_mu(list(mu = mu))
    expect_equal(p_mu, p_mu_ref)

    f_pi <- function(x) { dexp(x$pi, rate = 12.123, log = TRUE) }
    p_pi <- cpp_prior_pi(param, config, f_pi)
    p_pi_ref <- f_pi(list(pi = pi))
    expect_equal(p_pi, p_pi_ref)

})  
