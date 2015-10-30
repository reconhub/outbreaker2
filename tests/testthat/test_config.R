context("Test outbreaker data and settings")


## test data ##
test_that("test: data are processed fine", {
    ## skip on CRAN
    skip_on_cran()

    ## get data
    x <- fakeOutbreak
    out <- outbreaker.data(dates=x$collecDates, dna=x$dat$dna, w.dens=x$w)

    ## check output
    expect_is(out, "list")
    expect_is(out$D, "matrix")
    expect_equal(out$MAX.RANGE, 11)
    expect_error(outbreaker.data())
    expect_error(outbreaker.data(dates=1, w.dens=c(0,-1)))
    expect_error(outbreaker.data(dates=1, w.dens=c(0,1), f.dens=c(0,-1)))
})


## test settings ##
test_that("test: settings are processed fine", {
    ## skip on CRAN
    skip_on_cran()

    ## get data
    x <- fakeOutbreak
    dat <- outbreaker.data(dates=x$collecDates, dna=x$dat$dna, w.dens=x$w)

    ## check output
    expect_is(outbreaker.config(), "list")
    expect_is(outbreaker.config(data=dat), "list")
    expect_error(outbreaker.config(unknownArg))
    expect_error(outbreaker.config(init.tree="wrongtreeinit"))
    expect_error(outbreaker.config(init.mu=-5))
    expect_error(outbreaker.config(n.iter=0))
    expect_error(outbreaker.config(sample.every=0))
})


