context("Test prior functions")



## test outbreaker_create_priors ##
test_that("create_priors gives expected results", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    config <- create_config()
    out <- create_priors(config)

    ## tests
    expect_is(out, "list")
    expect_equal(length(out), 3)
    expect_equal(names(out), c("mu",
                               "pi",
                               "all")
                 )
    ## check that all items are functions
    expect_true(all(vapply(out, is.function, logical(1))))

    ## check that closure worked
    expect_identical(config$prior_mu, environment(out$mu)$rate)
    expect_identical(config$prior_pi,
                     c(environment(out$pi)$shape1, environment(out$pi)$shape2))

})





## test plausible values ##
test_that("priors have expected values", {
    ## skip on CRAN
    skip_on_cran()


    ## generate inputs
    config <- create_config()
    param <- outbreaker_create_mcmc(config,
                                    data = outbreaker_data())$current

    priors <- create_priors(config)

    ## generate outputs
    out_mu <- priors$mu(param)
    out_pi <- priors$pi(param)
    out_all <- priors$all(param)

    ## tests
    expect_equal(out_mu + out_pi, out_all)
    expect_true(!any(is.na(out_mu)))
    expect_true(!any(is.na(out_pi)))
    expect_true(!any(is.na(out_all)))
    expect_equal(out_all, 8.16209573106)

})






test_that("Prior customisation", {

    ## skip on CRAN
    skip_on_cran()

    
    ## check default config
    config <- create_config()
    expect_equal(create_priors(),
                     create_priors(config))

    ## check errors
    msg <- "The following priors are not functions: mu"
    expect_error(create_priors(mu = "Chtulhu"), msg)

    msg <- "The following priors dont' have a single argument: mu"
    expect_error(create_priors(mu = plot), msg)


    ## custom prior parameters
    config <- create_config(prior_pi = c(2, 1))
    p1 <- create_priors(config) # passing config
    p2 <- create_priors(config = list(prior_pi = c(1,1))) # passing list
    expect_equal(p1, p2)
    expect_equal(p1$pi(list(pi = 0.001)),
                 dbeta(0.001, 2,1, log = TRUE))
    
    ## custom function
    f_mu <- function(x){dexp(x$mu, rate = 100, log = TRUE)}
    expect_equal(create_priors(mu = f_mu)$mu, f_mu)
    
})  
