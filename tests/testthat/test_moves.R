context("Test movement of parameters and augmented data")

## test various movements  ##
test_that("parameters and augmented data move", {
    ## skip on CRAN
    skip_on_cran()

    ## generate data
    data(fakeOutbreak)
    data <- with(fakeOutbreak, outbreaker.data(dates=collecDates, w.dens=w, dna=dat$dna))
    config <- outbreaker.config(data=data)
    chain <- outbreaker.mcmc.init(data=data, config=config)
    rand <- outbreaker.rand.vec(config=config)

    ## ## test mu ##
    ## expect_equal(config$init.mu, move.mu(chain=chain, data=data, r.acc=0, r.new=1))
    ## expect_equal(config$init.mu+0.002, move.mu(chain=chain, data=data, r.acc=-Inf, r.new=0.002))

    ## test move.swap.ances ##
    res <- move.swap.ances(data=data, chain=chain, config=config)
    expect_equal(length(chain), length(res))
    expect_equal(length(unlist(chain)), length(unlist(res)))
    expect_equal(names(chain), names(res))
    res <- move.swap.ances(data=data, chain=chain, config=config)
    expect_equal(chain,res)
})

