context("Test outbreaker config")


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
    expect_equal(length(outbreaker.config(init.tree="random",data=dat)$ances), dat$N)
    expect_error(outbreaker.config(uknownarg=123))
    expect_error(outbreaker.config(init.tree="wrongtreeinit"))
    expect_error(outbreaker.config(init.mu=-5))
    expect_error(outbreaker.config(n.iter=0))
    expect_error(outbreaker.config(sample.every=0))
    expect_error(outbreaker.config(init.tree=1:5, data=dat))
    expect_warning(outbreaker.config(init.tree="seqTrack", data=dat.nodna))
    expect_warning(outbreaker.config(init.tree=rep(-1,dat$N), data=dat))
})


