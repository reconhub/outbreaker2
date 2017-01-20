context("Test detection of imported cases")

## test various movements  ##
test_that("Test detection of imported cases", {
    ## skip on CRAN
    skip_on_cran()


    ## generate inputs
    data(fake_outbreak)
    data <- with(fake_outbreak, outbreaker_data(dates = collecDates, w_dens = w, dna = dat$dna))
    config <- create_config(data = data)
    data <- add_convolutions(data = data, config = config)
    temp <- outbreaker_create_mcmc(data = data, config = config)
    param_current <- temp$current
    param_store <- temp$store

    ll <- create_loglike()
    priors <- create_priors(config)
    post <- create_posteriors(ll, priors)
    densities <- list(loglike = ll, priors = priors, posteriors = post)

    moves <- create_moves(config = config, data = data, densities = densities)

    ## detect imported cases
    out <- outbreaker_find_imports(moves = moves, data = data,
                                   param_current = param_current,
                                   param_store = param_store,
                                   config = config, densities = densities)

    ## tests ##
    expect_identical(which(!out$config$move_alpha), which(!out$config$move_kappa))
    expect_identical(out$param_store$alpha[[1]], out$param_current$alpha)
    expect_identical(out$param_store$kappa[[1]], out$param_current$kappa)
    expect_equal(which(is.na(out$param_current$alpha)), c(1,4,28))
    expect_true(all(config$move_alpha==!is.na(param_current$alpha)))

})


