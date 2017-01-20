context("Test movements")

## test various movements  ##
test_that("parameters and augmented data move", {
    ## skip on CRAN
    skip_on_cran()


    ## generate inputs
    data(fake_outbreak)
    data <- with(fake_outbreak, outbreaker_data(
        dates = collecDates, w_dens = w, dna = dat$dna))
    config <- create_config(data = data)
    
    config_no_move <- create_config(move_alpha = FALSE,
                                    move_t_inf = FALSE,
                                    move_mu = FALSE, move_pi = FALSE,
                                    move_kappa = FALSE, data = data)

    data <- add_convolutions(data = data, config = config)
    param <- outbreaker_create_mcmc(data = data, config = config)$current
    ll <- custom_likelihood()
    priors <- custom_priors(config)
    densities <- list(loglike = ll, priors = priors)

    moves <- custom_moves(config = config, data = data, densities = densities)
    moves_no_move <- custom_moves(config = config_no_move, densities = densities)


    ## test moves lists ##
    expect_equal(length(moves_no_move), 0L)
    expect_equal(length(moves), 6L)
    expect_true(all(vapply(moves, is.function, logical(1))))

    ## test moves ##
    for (i in seq_along(moves)) {
        ## test closures
        expect_identical(environment(moves[[i]])$config, config)
        expect_identical(environment(moves[[i]])$densities, densities)

        ## make moves
        set.seed(1)
        res <- moves[[i]](param = param)

        ## check that content in param after movements has identical shape
        expect_equal(length(param), length(res))
        expect_equal(length(unlist(param)), length(unlist(res)))
        expect_equal(names(param), names(res))

    }

})


