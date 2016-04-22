context("Test prior functions")



## test outbreaker.create.priors ##
test_that("create.priors gives expected results", {
    ## skip on CRAN
    skip_on_cran()
    

    ## generate data
    config <- outbreaker.config()
    out <- create.priors(config)

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
    expect_identical(config$prior.mu, environment(out$mu)$rate)
    expect_identical(config$prior.pi,
                     c(environment(out$pi)$shape1, environment(out$pi)$shape2))

})





## test plausible values ##
test_that("priors have expected values", {
    ## skip on CRAN
    skip_on_cran()
    

    ## generate inputs
    config <- outbreaker.config()
    param <- outbreaker.create.mcmc(config,
                                    data=outbreaker.data())

    priors <- create.priors(config)

    ## generate outputs
    out.mu <- priors$mu(param)
    out.pi <- priors$pi(param)
    out.all <- priors$all(param)

    ## tests
    expect_equal(out.mu + out.pi, out.all)
    expect_true(!any(is.na(out.mu)))
    expect_true(!any(is.na(out.pi)))
    expect_true(!any(is.na(out.all)))
    expect_equal(out.all, 8.16209573106)

})
