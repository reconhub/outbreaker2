context("Test binding functions")

## test output format ##
test_that("", {
    ## skip on CRAN
    skip_on_cran()


    ## get data
    data(fake_outbreak)
    dat <- fake_outbreak$dat
    w <- fake_outbreak$w

    ## run outbreaker
    data <- list(dna = dat$dna, dates = dat$onset, w_dens = w)
    config <- list(n_iter = 10, sample_every = 1, paranoid = TRUE,
                   find_import = FALSE)

})
