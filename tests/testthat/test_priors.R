context("Test prior functions")


## test plausible values ##
test_that("priors have plausible values", {
    ## skip on CRAN
    skip_on_cran()
    rm(list=ls())

    ## generate data
    param <- outbreaker.create.mcmc(config=outbreaker.config(),
                                    data=outbreaker.data())


    ## tests
    out.mu <- priors$mu(param)
    out.pi <- priors$pi(param)
    out.all <- priors$all(param)
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
