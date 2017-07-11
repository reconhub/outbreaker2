context("Test sim_ctd")

## test output format ##
test_that("sim_ctd working as expected", {
    ## skip on CRAN
    skip_on_cran()

    ## get data
    data(fake_outbreak)
    x <- fake_outbreak
    N <- length(x$ances)
    tTree <- data.frame(i = x$ances, j = 1:N)

    ## error check
    expect_error(sim_ctd(tTree = tTree, eps = 1.1, lambda = 0),
                 "eps and lambda must be probabilities")
    
    ## for eps = 1 / lambda = 0, there should be one
    ## contact per transmission pair
    ctd <- sim_ctd(tTree = tTree, eps = 1, lambda = 0)
    expect_equal(nrow(ctd), sum(!is.na(x$ances)))

    ## for eps = 1 / lambda = 1, ensure all possible
    ## contacts are reported
    ctd <- sim_ctd(tTree = tTree, eps = 1, lambda = 1)
    expect_equal(nrow(ctd), N*(N-1)/2)

})
