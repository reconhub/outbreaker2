context("Test prior functions")


## test plausible values ##
test_that("priors have plausible values", {
    ## skip on CRAN
    skip_on_cran()
    rm(list=ls())

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

    expect_true(!any(is.na(out)))

})



## test outbreaker.create.priors ##
test_that("outbreaker.create.priors gives expected results", {
    ## skip on CRAN
    skip_on_cran()
    rm(list=ls())

    ## generate data
    out <- outbreaker.create.priors()

    ## tests
    expect_is(out, "list")
    expect_equal(length(out), 2)
    expect_equal(names(out), c("mu",
                               "all")
                 )

})
