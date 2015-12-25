context("Test outbreaker data and settings")


## test data ##
test_that("test: data are processed fine", {
    ## skip on CRAN
    skip_on_cran()

    ## get data
    x <- fakeOutbreak
    out <- outbreaker.data(dates=x$collecDates, dna=x$dat$dna, w.dens=x$w)
    out.nodna <- outbreaker.data(dates=x$collecDates, w.dens=x$w)

    ## check output
    expect_is(out, "list")
    expect_is(out$D, "matrix")
    expect_equal(out$MAX.RANGE, 11)
    expect_equal(out.nodna$L, 0)
    expect_equal(out$L, 1e4)
    expect_error(outbreaker.data(dates=1, w.dens=c(0,-1)))
    expect_error(outbreaker.data(dates=1, w.dens=c(0,1), f.dens=c(0,-1)))
})

