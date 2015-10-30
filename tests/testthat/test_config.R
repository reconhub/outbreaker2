context("Test outbreaker settings")

## test output format ##
test_that("test: settings are processed fine", {
    ## skip on CRAN
    skip_on_cran()

    ## check output
    expect_is(outbreaker.config(), "list")
    expect_error(outbreaker.config(unknownArg))
    expect_error(outbreaker.config(init.tree="wrongtreeinit"))
    expect_error(outbreaker.config(init.mu=-5))
    expect_error(outbreaker.config(n.iter=0))
    expect_error(outbreaker.config(sample.every=0))
})

