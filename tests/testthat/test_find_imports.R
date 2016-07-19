context("Test detection of imported cases")

## test various movements  ##
test_that("Test detection of imported cases", {
    ## skip on CRAN
    skip_on_cran()


    ## generate inputs
    data(fake.outbreak)
    data <- with(fake.outbreak, outbreaker.data(dates=collecDates, w.dens=w, dna=dat$dna))
    config <- outbreaker.config(data=data)
    data <- add.convolutions(data=data, config=config)
    temp <- outbreaker.create.mcmc(data=data, config=config)
    param.current <- temp$current
    param.store <- temp$store

    ll <- create.loglike(data)
    priors <- create.priors(config)
    post <- create.posteriors(ll, priors)
    densities <- list(loglike=ll, priors=priors, posteriors=post)

    moves <- create.moves(config=config, densities=densities)

    ## detect imported cases
    out <- outbreaker.find.imports(moves = moves, data = data,
                                   param.current = param.current,
                                   param.store = param.store,
                                   config = config, densities = densities)

    ## tests ##
    expect_identical(which(!out$config$move.alpha), which(!out$config$move.kappa))
    expect_identical(out$param.store$alpha[[1]], out$param.current$alpha)
    expect_identical(out$param.store$kappa[[1]], out$param.current$kappa)
    expect_equal(which(is.na(out$param.current$alpha)), c(1,4,28))
    expect_true(all(config$move.alpha==!is.na(param.current$alpha)))

})


