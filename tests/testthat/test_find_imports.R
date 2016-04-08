context("Test detection of imported cases")

## test various movements  ##
test_that("Test detection of imported cases", {
    ## skip on CRAN
    skip_on_cran()
    rm(list=ls())

    ## generate data
    data(fakeOutbreak)
    data <- with(fakeOutbreak, outbreaker.data(dates=collecDates, w.dens=w, dna=dat$dna))
    config <- outbreaker.config(data=data)
    param <- outbreaker.mcmc.init(data=data, config=config)
    rand <- outbreaker.rand.vec(config=config)
    moves <- outbreaker.create.moves(config=config)
    data <- add.convolutions(data=data, config=config)

    ## detect imported cases
    out <- outbreaker.find.imports(moves=moves, data=data, config=config, param=param, rand=rand)

    ## tests ##
    expect_identical(which(!out$config$move.alpha), which(!out$config$move.kappa))
    expect_identicalidentical(out$param$alpha[[1]], out$param$current.alpha)
    expect_identical(out$param$kappa[[1]], out$param$current.kappa)
    expect_equal(which(is.na(out$param$current.alpha)), c(1,4,28))
    expect_true(all(config$move.alpha==!is.na(param$current.alpha)))
})


