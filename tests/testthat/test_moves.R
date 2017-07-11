context("Test movements")

## test various movements  ##
test_that("Movements preserve param structure", {
    ## skip on CRAN
    skip_on_cran()


    ## generate inputs
    data(fake_outbreak)
    data <- with(fake_outbreak,
                 outbreaker_data(dates = onset,
                                 w_dens = w,
                                 dna = dna))
    config <- create_config(data = data)

    config_no_move <- create_config(
      move_alpha = FALSE,
      move_t_inf = FALSE,
      move_mu = FALSE, move_pi = FALSE,
      move_eps = FALSE, move_lambda = FALSE,
      move_kappa = FALSE,
      move_swap_cases = FALSE, data = data)

    data <- add_convolutions(data = data, config = config)
    param <- create_param(data = data, config = config)$current
    ll <- custom_likelihoods()
    priors <- custom_priors()

    moves <- bind_moves(config = config, data = data,
                          likelihoods = ll, priors = priors)
    moves_no_move <- bind_moves(config = config_no_move,
                                  likelihoods = ll, priors = priors)


    ## test moves lists ##
    expect_equal(length(moves_no_move), 0L)
    expect_equal(length(moves), 6L)
    expect_true(all(vapply(moves, is.function, logical(1))))




    ## test moves ##
    for (i in seq_along(moves)) {

        ## chech closure: data
        expect_identical(environment(moves[[i]])$data, data)

        ## make moves
        set.seed(1)
        res <- moves[[i]](param = param)

        ## check that content in param after movements has identical shape
        expect_equal(length(param), length(res))
        expect_equal(length(unlist(param)), length(unlist(res)))
        expect_equal(names(param), names(res))

    }

})






test_that("Binding of moves works", {
    ## skip on CRAN
    skip_on_cran()


    ## generate inputs
    data(fake_outbreak)
    data <- with(fake_outbreak,
                 outbreaker_data(dates = onset,
                                 w_dens = w,
                                 dna = dna))
    config <- create_config(data = data)
    data <- add_convolutions(data = data, config = config)
    param <- create_param(data = data, config = config)$current
    ll <- custom_likelihoods()
    priors <- custom_priors()

    ## check re-input consistency
    expect_identical(custom_moves(),
                     custom_moves(custom_moves()))

    ## check custom_moves defaults
    moves <- custom_moves()

    expect_length(moves, 8L)
    expect_true(all(vapply(moves, is.function, FALSE)))
    expect_named(moves)
    expected_names <- c("mu", "pi", "eps", "lambda", "alpha", "swap_cases", "t_inf", "kappa")
    expect_true(all(expected_names %in% names(moves)))


    ## check binding
    moves <- bind_moves(moves, config = config, data = data,
                          likelihoods = ll, priors = priors)

    exp_names <- c("custom_prior", "custom_ll", "config", "data")
    expect_true(all(exp_names %in% names(environment(moves$mu))))

    exp_names <- c("custom_prior", "custom_ll", "config", "data")
    expect_true(all(exp_names %in% names(environment(moves$pi))))

    exp_names <- c("list_custom_ll", "data")
    expect_true(all(exp_names %in% names(environment(moves$alpha))))

    exp_names <- c("list_custom_ll", "data")
    expect_true(all(exp_names %in% names(environment(moves$swap_cases))))

    exp_names <- c("list_custom_ll", "data")
    expect_true(all(exp_names %in% names(environment(moves$t_inf))))

    exp_names <- c("list_custom_ll", "config", "data")
    expect_true(all(exp_names %in% names(environment(moves$kappa))))

})






test_that("Customisation of moves works", {
    ## skip on CRAN
    skip_on_cran()


    ## generate inputs
    data(fake_outbreak)
    data <- with(fake_outbreak,
                 outbreaker_data(dates = onset,
                                 w_dens = w,
                                 dna = dna))
    config <- create_config(data = data, n_iter = 1000,
                            find_import = FALSE,
                            sample_every = 10)
    data <- add_convolutions(data = data, config = config)
    param <- create_param(data = data, config = config)$current
    ll <- custom_likelihoods()
    priors <- custom_priors()


    ## check custom movement for mu - outside outbreaker
    f <- function(param, data, config = NULL) {
        return(param)
    }

    moves <- bind_moves(list(mu = f), config = config, data = data,
                          likelihoods = ll, priors = priors)

    expect_identical(body(moves$mu), body(f))
    expect_identical(names(formals(moves$mu)), "param")
    expect_identical(data, environment(moves$mu)$data)
    expect_identical(config, environment(moves$mu)$config)
    expect_identical(moves$mu(param), param)


    ## same check, run within outbreaker
    out <- outbreaker(data, config, moves = list(mu = f))
    expect_true(all(out$mu == 1e-4))

})






## test swapping and temporal ordering  ##
test_that("Swap equally likely index cases", {
    ## skip on CRAN
    skip_on_cran()


    ## generate inputs
    data <- outbreaker_data(dates = c(50, 51, 110),
                            w_dens = rep(1, 100))
    config <- create_config(init_kappa = 1,
                            move_kappa = FALSE,
                            find_import = FALSE,
                            data = data)

    set.seed(1)
    res <- outbreaker(data, config)
    table(res$alpha_1)
    table(res$alpha_2)
    table(res$alpha_3)

})
