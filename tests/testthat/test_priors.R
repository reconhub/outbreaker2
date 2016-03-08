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
