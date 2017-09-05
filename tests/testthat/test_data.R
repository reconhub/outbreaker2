context("Test outbreaker data and settings")


## test data ##
test_that("test: data are processed fine", {
    ## skip on CRAN
    skip_on_cran()


    ## get data
    x <- fake_outbreak
    out <- outbreaker_data(dates = x$onset, dna = x$dna, w_dens = x$w)
    out_nodna <- outbreaker_data(dates = x$onset, w_dens = x$w)

    ## check output
    expect_is(out, "list")
    expect_is(out$D, "matrix")
    expect_equal(out$max_range, 11)
    expect_equal(out_nodna$L, 0)
    expect_equal(out$L, 1e4)
    expect_equal(out$w_dens, out$f_dens)
    expect_equal(out$log_w_dens[1,], out$log_f_dens)
    expect_warning(outbreaker_data(dates = 1, w_dens = c(0, 0.5, 0.1, 0, 0)),
                   "Removed trailing zeroes found in w_dens")
    expect_error(outbreaker_data(dates = 1, w_dens = c(0,-1)),
                 "w_dens has negative entries")

    expect_error(outbreaker_data(dates = 1, w_dens = c(0,1), f_dens = c(0,-1)),
                 "f_dens has negative entries")

    wrong_lab_dna <- x$dna
    rownames(wrong_lab_dna) <- paste0("host_", seq_len(nrow(wrong_lab_dna)))
    expect_error(outbreaker_data(dates = x$onset, dna = wrong_lab_dna, w_dens = x$w),
                 "DNA sequence labels don't match case ids")


})

