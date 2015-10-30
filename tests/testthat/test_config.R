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
    dat.nodna <- dat
    dat.nodna$dna <- NULL

    ## check output
    expect_is(outbreaker.config(), "list")
    expect_is(outbreaker.config(data=dat), "list")
    expect_equal(outbreaker.config(init.tree="star", data=dat)$ances, c(NA, rep(1,29)))
    expect_error(outbreaker.config(uknownarg=123))
    expect_error(outbreaker.config(init.tree="wrongtreeinit"))
    expect_error(outbreaker.config(init.mu=-5))
    expect_error(outbreaker.config(n.iter=0))
    expect_error(outbreaker.config(sample.every=0))
    expect_error(outbreaker.config(init.tree=1:5, data=dat))
    expect_warning(outbreaker.config(init.tree="seqTrack", data=dat.nodna))
    expect_warning(outbreaker.config(init.tree=rep(-1,dat$N), data=dat))
})


