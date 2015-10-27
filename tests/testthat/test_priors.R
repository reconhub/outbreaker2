context("Test prior functions")


## test plausible values ##
test_that("priors have plausible values", {
    ## skip on CRAN
    skip_on_cran()

    ## generate data
    mu <- runif(100)

    ## tests
    prior.mu <- sapply(mu, prior.mu)
    expect_true(!any(is.na(prior.mu)))
})
