context("Test outbreaker data and settings")


## test data ##
test_that("test: data are processed fine", {
    ## skip on CRAN
    skip_on_cran()


    ## get data
    x <- fake_outbreak
    out <- outbreaker_data(dates = x$collecDates, dna = x$dat$dna, w_dens = x$w)
    out_nodna <- outbreaker_data(dates = x$collecDates, w_dens = x$w)

    ## check output
    expect_is(out, "list")
    expect_is(out$D, "matrix")
    expect_equal(out$max_range, 11)
    expect_equal(out_nodna$L, 0)
    expect_equal(out$L, 1e4)
    expect_equal(out$w_dens, out$f_dens)
    expect_equal(out$log_w.dens[1,], out$log_f.dens)
    expect_error(outbreaker_data(dates = 1, w_dens = c(0,-1)),
                 "w_dens has negative entries")

    expect_error(outbreaker_data(dates = 1, w_dens = c(0,1), f_dens = c(0,-1)),
                 "f_dens has negative entries")

})

