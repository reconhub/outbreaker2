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
                      config=list(n.iter=100, sample.every=10, safemode=TRUE))
    out.df <- as.data.frame(out)

    ## check output
    expect_is(out, "mcmc")
    expect_is(out.df, "data.frame")
    expect_equal(nrow(out), 11)
    expect_true(!any(is.na(out.df$post)))
    expect_true(all(out.df$post> -1e30))
})

