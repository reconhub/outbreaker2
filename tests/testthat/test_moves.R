context("Test movements")

## test various movements  ##
test_that("parameters and augmented data move", {
    ## skip on CRAN
    skip_on_cran()
    

    ## generate inputs
    data(fakeOutbreak)
    data <- with(fakeOutbreak, outbreaker.data(dates=collecDates, w.dens=w, dna=dat$dna))
    config <- outbreaker.config(data=data)
    config.no.move <- outbreaker.config(data=data, move.alpha=FALSE, move.t.inf=FALSE,
                                        move.mu=FALSE, move.pi=FALSE, move.kappa=FALSE)
    data <- add.convolutions(data=data, config=config)
    param <- outbreaker.create.mcmc(data=data, config=config)
    ll <- create.loglike(data)
    priors <- create.priors(config)
    post <- create.posteriors(ll, priors)
    densities <- list(loglike=ll, priors=priors, posteriors=post)

    rand <- outbreaker.rand.vec(config=config)
    moves <- create.moves(config=config, densities=densities, rand=rand)
    moves.no.move <- create.moves(config=config.no.move, densities=densities, rand=rand)


    ## test moves lists ##
    expect_equal(length(moves.no.move), 0L)
    expect_equal(length(moves), 6L)
    expect_true(all(vapply(moves, is.function, logical(1))))

    ## test moves ##
    for (i in seq.int(length(moves))) {
        ## test closures
        identical(environment(moves[[i]])$config, config)
        identical(environment(moves[[i]])$densities, densities)
        identical(environment(moves[[i]])$rand, rand)

        ## make moves
        set.seed(1)
        res <- moves[[i]](param=param)

        ## check that content in param after movements has identical shape
        expect_equal(length(param), length(res))
        expect_equal(length(unlist(param)), length(unlist(res)))
        expect_equal(names(param), names(res))

    }

})


