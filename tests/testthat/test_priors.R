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
