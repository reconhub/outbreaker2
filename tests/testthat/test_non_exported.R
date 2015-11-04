context("Test non-exported functions")


## test data ##
test_that("test: find.possible.ances", {
    ## skip on CRAN
    skip_on_cran()

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

