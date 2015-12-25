context("Test prior functions")


## test plausible values ##
test_that("priors have plausible values", {
    ## skip on CRAN
    skip_on_cran()
    rm(list=ls())

    ## generate data
    param <- lapply(1:30, function(i) outbreaker.mcmc.init(
        config=outbreaker.config(init.mu=runif(1)),
        data=outbreaker.data()))

    ## tests
    out <- sapply(param, prior.mu)
    expect_true(!any(is.na(out)))
})
