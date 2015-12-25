context("Test movement of parameters and augmented data")

## test various movements  ##
test_that("parameters and augmented data move", {
    ## skip on CRAN
    skip_on_cran()
    rm(list=ls())

    ## generate data
    data(fakeOutbreak)
    data <- with(fakeOutbreak, outbreaker.data(dates=collecDates, w.dens=w, dna=dat$dna))
    config <- outbreaker.config(data=data)
    param <- outbreaker.mcmc.init(data=data, config=config)
    rand <- outbreaker.rand.vec(config=config)

    ## test move.swap.ances ##
    res <- move.swap.ances(data=data, param=param, config=config, rand=rand)
    expect_equal(length(param), length(res))
    expect_equal(length(unlist(param)), length(unlist(res)))
    expect_equal(names(param), names(res))
})

