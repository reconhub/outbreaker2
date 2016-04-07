context("Test outbreaker")

## test output format ##
test_that("test: outbreaker's output have expected format", {
    ## skip on CRAN
    skip_on_cran()
    rm(list=ls())

    ## get data
    data(fakeOutbreak)
    dat <- fakeOutbreak$dat
    w <- fakeOutbreak$w

    ## run outbreaker
    out <- outbreaker(data=list(dna=dat$dna, dates=dat$onset, w.dens=w),
                      config=list(n.iter=100, sample.every=10, paranoid=TRUE))
    out.df <- as.data.frame(out)

    out2 <- outbreaker(data=list(dates=dat$onset, w.dens=w),
                      config=list(n.iter=100, sample.every=10, paranoid=TRUE))

    ## check output
    expect_is(out, "outbreaker.chains")
    expect_is(out.df, "data.frame")
    expect_equal(nrow(out), 11)
    expect_true(!any(is.na(out.df$post)))
    expect_true(all(out.df$post> -1e30))

})




## test convergence ##
test_that("test: convergence to decent results for toy example", {
    ## skip on CRAN
    skip_on_cran()
    rm(list=ls())

    ## get data
    data(fakeOutbreak)
    dat <- fakeOutbreak$dat
    w <- fakeOutbreak$w

    ## outbreaker DNA + time
    set.seed(1)
    out <- outbreaker(data=list(dna=dat$dna, dates=dat$onset, w.dens=w),
                      config=list(n.iter=5000, sample.every=100, init.tree="star"))

    ## check that tree is OK
    out.smry <- summary(out, burn=1000)
    expect_true(min(out.smry$post) > -920) # approx log post values
    expect_true(mean(out.smry$tree$from==dat$ances, na.rm=TRUE) > .85) # at least 85% ancestries correct
    expect_true(mean(abs(out.smry$tree$time - dat$onset), na.rm=TRUE)<3) # infection datewithin 3 days on average
    expect_true(min(out.smry$mu) > 0.0002 && max(out.smry$mu) < 0.0004) # mu between 2e-4 and 4 e-4

    ## outbreaker time, no DNA
    set.seed(1)
    out.no.dna <- outbreaker(data=list(dates=dat$onset, w.dens=w),
                      config=list(n.iter=5000, sample.every=100, init.tree="star"))

    ## check that tree is OK
    out.no.dna.smry <- summary(out.no.dna, burn=1000)
    expect_true(min(out.no.dna.smry$post) > -100) # approx log post values
    expect_true(mean(out.no.dna.smry$tree$from==dat$ances, na.rm=TRUE) > .05) # at least 5% ancestries correct
    expect_true(mean(abs(out.no.dna.smry$tree$time - dat$onset), na.rm=TRUE)<3) # infection datewithin 3 days on average

})
