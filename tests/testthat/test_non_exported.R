context("Test non-exported functions")


## test find.possible.ances ##
test_that("test: find.possible.ances", {
    ## skip on CRAN
    skip_on_cran()
    rm(list=ls())

    ## get data
    ans1 <- outbreaker2:::find.possible.ances(1:10, 1)
    ans2 <- outbreaker2:::find.possible.ances(1:10, 10)
    ans3 <- outbreaker2:::find.possible.ances(1:10, 2)

    ## check output
    expect_true(is.na(ans1))
    expect_is(ans2, "integer")
    expect_true(ans2<10 || ans2>0)
    expect_equal(ans3, 1L)
    expect_error(find.possible.ances(1:10, NA))
})


## test swap.ances ##
test_that("Auxiliary functions for ancestries are working", {
    ## skip on CRAN
    skip_on_cran()

    ## generate data
    ances <- c(NA,1,1,2,2)
    t.inf <- c(0, 2,3, 5,7)

    ## tests non swapping (imported case)
    expect_warning(res <- outbreaker2:::swap.ances(ances, t.inf, 1))
    expect_equal(ances, res$ances)
    expect_equal(t.inf, res$t.inf)

    ## test full swapping
    res <- outbreaker2:::swap.ances(ances, t.inf, 2)
    expect_equal(res$ances, c(2,NA,2,1,1))
    expect_equal(res$t.inf, c(2,0,3,5,7))

})
