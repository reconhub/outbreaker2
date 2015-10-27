context("Test outbreaker")

data(fakeOutbreak)
dat <- fakeOutbreak$dat
w <- fakeOutbreak$w

## test output format ##
test_that("test: outbreaker's output have expected format", {
    ## skip on CRAN
    skip_on_cran()

    ## run outbreaker
    out <- outbreaker(dna=dat$dna, dates=dat$onset, w.dens=w, n.iter=10)
    out.df <- as.data.frame(out)

    ## check output
    expect_is(out, "mcmc")
    expect_is(out.df, "data.frame")
    expect_equal(nrow(out), 10)
    expect_true(!any(is.na(out.df$post)))
    expect_true(all(out.df$post> -1e30))
})

