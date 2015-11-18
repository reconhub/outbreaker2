context("Test movement of parameters and augmented data")


## test plausible values ##
test_that("Auxiliary functions for ancestries are working", {
    ## skip on CRAN
    skip_on_cran()

    ## generate data
    ances <- c(NA,1,1,2,2)
    t.inf <- c(0, 2,3, 5,7)

    ## test swap.ances
    expect_warning(res <- swap.ances(ances, t.inf, 1))
    expect_equal(ances, res$ances)
    expect_equal(t.inf, res$t.inf)

})
