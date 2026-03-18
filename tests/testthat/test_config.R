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
    expect_error(create_config(pb = 'TRUE'),
                 "pb must be a logical")


})




test_that("validation does not alter valide objects", {
    c1 <- create_config()
    c2 <- create_config(c1)
    expect_identical(c1, c2)
})

test_that("move_alpha recycles only when length 1, else must match N", {
    skip_on_cran()
    x <- fake_outbreak
    dat <- outbreaker_data(dates = x$sample, dna = x$dna, w_dens = x$w)
    ## single value recycles to N
    config1 <- create_config(data = dat, move_alpha = TRUE)
    expect_length(config1$move_alpha, dat$N)
    expect_true(all(config1$move_alpha))
    config1b <- create_config(data = dat, move_alpha = FALSE)
    expect_length(config1b$move_alpha, dat$N)
    expect_false(any(config1b$move_alpha))
    ## full vector of length N is used as-is (no recycle); imports still set FALSE
    move_alpha_N <- rep(c(TRUE, FALSE), length.out = dat$N)
    config2 <- create_config(data = dat, move_alpha = move_alpha_N)
    expect_length(config2$move_alpha, dat$N)
    ## non-imports keep the supplied value; imports are forced FALSE later in create_config
    non_import <- !is.na(config2$init_alpha)
    expect_equal(config2$move_alpha[non_import], move_alpha_N[non_import])
    ## wrong length errors
    expect_error(create_config(data = dat, move_alpha = c(TRUE, FALSE)),
                 "move_alpha must be of length 1 or N")
})
