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
    alpha <- rep(1, x$dat$n)

    ## check output
    expect_is(outbreaker.config(), "list")
    expect_is(outbreaker.config(data=dat), "list")
    expect_equal(outbreaker.config(init.tree="star", data=dat)$init.alpha, c(NA, rep(1,29)))
    expect_equal(outbreaker.config(init.tree=alpha)$init.tree, outbreaker.config(init.tree=alpha)$init.alpha)
    expect_equal(length(outbreaker.config(init.tree="random",data=dat)$init.alpha), dat$N)
    expect_equal(sort(unique(outbreaker.config(data=dat, init.kappa=1:2)$init.kappa)), 1:2)
    expect_error(outbreaker.config(uknownarg=123),
                 "Additional invalid options: uknownarg")
    expect_error(outbreaker.config(init.tree="wrongtreeinit"),
                 "'arg' should be one of “seqTrack”, “star”, “random”")
    expect_error(outbreaker.config(init.mu=-5),
                 "init.mu is negative")
    expect_error(outbreaker.config(n.iter=0),
                 "n.iter is smaller than 2")
    expect_error(outbreaker.config(sample.every=0),
                 "sample.every is smaller than 1")
    expect_error(outbreaker.config(init.tree=1:5, data=dat),
                 "inconvenient length for init.alpha")
    expect_message(outbreaker.config(init.tree="seqTrack", data=dat.nodna),
                   "Can't use seqTrack initialization with missing DNA sequences; using a star-like tree")
    expect_warning(outbreaker.config(init.tree=rep(-1,dat$N), data=dat),
                   "some initial ancestries refer to unknown cases")
})


