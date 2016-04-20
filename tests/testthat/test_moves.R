context("Test movement of parameters and augmented data")

## test various movements  ##
test_that("parameters and augmented data move", {
    ## skip on CRAN
    skip_on_cran()
    rm(list=ls())

    ## generate inputs
    data(fakeOutbreak)
    data <- with(fakeOutbreak, outbreaker.data(dates=collecDates, w.dens=w, dna=dat$dna))
    config <- outbreaker.config(data=data)
    config.no.move <- outbreaker.config(data=data, move.alpha=FALSE, move.t.inf=FALSE,
                                        move.mu=FALSE, move.pi=FALSE, move.kappa=FALSE)
    param <- outbreaker.create.mcmc(data=data, config=config)
    ll <- create.loglike(data)
    priors <- create.priors(config)
    post <- create.posteriors(ll, priors)
    densities <- list(loglike=ll, priors=priors, posteriors=post)

    rand <- outbreaker.rand.vec(config=config)
    moves <- outbreaker.create.moves(config=config)
    moves.no.move <- outbreaker.create.moves(config=config.no.move)

    ## test move.swap.cases ##
    res <- moves$move.swap.cases(param=param, config=config, densities=densities, rand=rand)

    expect_equal(length(param), length(res))
    expect_equal(length(unlist(param)), length(unlist(res)))
    expect_equal(names(param), names(res))
    expect_equal(length(moves.no.move), 0L)
})


