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
    out2 <- outbreaker(data=list(dates=dat$onset, w.dens=w),
                       config=list(n.iter=100, sample.every=10, paranoid=TRUE))

    ## outbreaker time only

    ## outbreaker DNA only



})
