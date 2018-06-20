context("Test outbreaker config")


## test settings ##
test_that("test: settings are processed fine", {
    ## skip on CRAN
    skip_on_cran()


    ## get data
    x <- fake_outbreak
    dat <- outbreaker_data(dates = x$sample, dna = x$dna, w_dens = x$w)
    dat_nodna <- dat
    dat_nodna$D <- dat_nodna$dna <- NULL
    alpha <- rep(1, length(x$ances))

    ## check output
    expect_is(create_config(), "list")
    expect_is(create_config(), "outbreaker_config")
    expect_is(create_config(data = dat), "list")
    expect_is(create_config(data = dat), "outbreaker_config")
    expect_equal(create_config(init_tree="star", data = dat)$init_alpha,
                 c(NA, rep(1,29)))
    expect_equal(create_config(init_tree = alpha)$init_tree,
                 create_config(init_tree = alpha)$init_alpha)
    expect_equal(length(create_config(init_tree="random",data = dat)$init_alpha), dat$N)
    expect_equal(sort(unique(create_config(data = dat, init_kappa = 1:2)$init_kappa)), 1:2)
    expect_error(create_config(uknownarg = 123),
                 "Additional invalid options: uknownarg")
    expect_error(create_config(init_tree="wrongtreeinit"),
                 "should be one of")
    expect_error(create_config(init_mu=-5),
                 "init_mu is negative")
    expect_error(create_config(n_iter = 0),
                 "n_iter is smaller than 2")
    expect_error(create_config(sample_every = 0),
                 "sample_every is smaller than 1")
    expect_error(create_config(init_tree = 1:5, data = dat),
                 "inconvenient length for init_alpha")
    expect_message(create_config(init_tree="seqTrack", data = dat_nodna),
                   "Can't use seqTrack initialization with missing DNA sequences; using a star-like tree")
    expect_warning(create_config(init_tree = rep(-1,dat$N), data = dat),
                   "some initial ancestries refer to unknown cases")
    expect_false(create_config(data = dat_nodna, move_mu = TRUE)$move_mu)

})




test_that("validation does not alter valide objects", {
    c1 <- create_config()
    c2 <- create_config(c1)
    expect_identical(c1, c2)
})
