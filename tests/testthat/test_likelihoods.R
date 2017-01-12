context("Test likelihood functions")



######################################################
######################################################
## End of the 'reference' functions. Tests are below.
## test ll_timing_infections ##
######################################################
######################################################

test_that("ll_timing_infections gives expected results", {
    ## skip on CRAN
    skip_on_cran()

    ## generate data
    times <- 0:4
    alpha <- c(NA,rep(1,4))
    w <- c(.1, .2, .5, .2, .1)
    data <- outbreaker_data(dates = times, w_dens = w)
    config <- outbreaker_config(data = data, init_tree = alpha)
    ll <- outbreaker2:::create_loglike(data)
    param <- outbreaker_create_mcmc(data = data, config = config)$current
    few_cases <- as.integer(c(1,3,4))
    rnd_cases <- sample(sample(seq_len(data$N), 3, replace = FALSE))


    ## tests
    out <- ll$timing_infections(param)
    out_few_cases <- ll$timing_infections(param, few_cases)
    out_rnd_cases <- ll$timing_infections(param, rnd_cases)
    ref <- .ll_timing_infections(data, param)
    ref_few_cases <- .ll_timing_infections(data, param, few_cases)
    ref_rnd_cases <- .ll_timing_infections(data, param, rnd_cases)

    expect_is(out, "numeric")
    expect_equal(out, -6.214608098)
    expect_equal(out_few_cases, -2.30258509299405)

    ## test against reference
    expect_equal(out, ref)
    expect_equal(out_few_cases, ref_few_cases)
    expect_equal(out_rnd_cases, ref_rnd_cases)
})




## test ll$timing_sampling ##
test_that("ll$timing_sampling gives expected results", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    times <- 0:4
    alpha <- c(NA,rep(1,4))
    samp_times <- times + c(1, 1, 2, 3, 4)
    f <- c(.1, .2, .5, .2, .1)
    data <- outbreaker_data(dates = samp_times, f_dens = f)
    ll <- outbreaker2:::create_loglike(data)
    config <- outbreaker_config(data = data, init_t_inf = times, init_tree = alpha)
    param <- outbreaker_create_mcmc(data = data, config = config)$current
    few_cases <- as.integer(c(1,3,4))
    rnd_cases <- sample(sample(seq_len(data$N), 3, replace = FALSE))


    ## tests
    out <- ll$timing_sampling(param)
    out_few_cases <- ll$timing_sampling(param, few_cases)
    out_rnd_cases <- ll$timing_sampling(param, rnd_cases)
    ref <- .ll_timing_sampling(data, param)
    ref_few_cases <- .ll_timing_sampling(data, param, few_cases)
    ref_rnd_cases <- .ll_timing_sampling(data, param, rnd_cases)

    expect_is(out, "numeric")
    expect_equal(out, -8.51719319142)
    expect_equal(out_few_cases, -4.60517018598809)

    ## test against reference
    expect_equal(out, ref)
    expect_equal(out_few_cases, ref_few_cases)
    expect_equal(out_rnd_cases, ref_rnd_cases)

})




## test ll$genetic ##
test_that("ll$genetic gives expected results", {
    ## skip on CRAN ##
    skip_on_cran()

    ## generate data ##
    data(fake_outbreak)
    data <- with(fake_outbreak, outbreaker_data(dates = collecDates, w_dens = w, dna = dat$dna))
    ll <- outbreaker2:::create_loglike(data)
    config <- outbreaker_config(data = data, init_mu = 0.543e-4)
    param <- outbreaker_create_mcmc(data = data, config = config)$current
    few_cases <- as.integer(c(1,3,4))
    rnd_cases <- sample(sample(seq_len(data$N), 5, replace = FALSE))

    ## tests ##
    ## expected values
    out <- ll$genetic(param)
    out_few_cases <- ll$genetic(param, few_cases)
    out_rnd_cases <- ll$genetic(param, rnd_cases)
    ref <- .ll_genetic(data, param)
    ref_few_cases <- .ll_genetic(data, param, few_cases)
    ref_rnd_cases <- .ll_genetic(data, param, rnd_cases)

    expect_is(out, "numeric")
    expect_equal(out, -997.840630501522)
    expect_equal(out_few_cases, -266.251194283819)
    

    ## test against R reference
    expect_equal(out, ref)
    expect_equal(out_few_cases, ref_few_cases)
    expect_equal(out_rnd_cases, ref_rnd_cases)

})






## test ll$reporting ##
test_that("ll$reporting gives expected results", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    data(fake_outbreak)
    data <- with(fake_outbreak, outbreaker_data(dates = collecDates, w_dens = w, dna = dat$dna))
    config <- outbreaker_config(data = data, init_mu = 0.543e-4)
    ll <- outbreaker2:::create_loglike(data)
    param <- outbreaker_create_mcmc(data = data, config = config)$current
    few_cases <- as.integer(c(1,3,4))
    rnd_cases <- sample(sample(seq_len(data$N), 3, replace = FALSE))


    ## tests
    out <- ll$reporting(param)
    out_few_cases <- ll$reporting(param, few_cases)
    out_rnd_cases <- ll$reporting(param, rnd_cases)
    ref <- .ll_reporting(data, param)
    ref_few_cases <- .ll_reporting(data, param, few_cases)
    ref_rnd_cases <- .ll_reporting(data, param, rnd_cases)

    expect_is(out, "numeric")
    expect_equal(out, -3.05545495407696)
    expect_equal(out_few_cases, -0.210721031315653)

    ## test against reference
    expect_equal(out, ref)
    expect_equal(out_few_cases, ref_few_cases)
    expect_equal(out_rnd_cases, ref_rnd_cases)

})






## test ll$timing ##
test_that("ll$timing gives expected results", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    data(fake_outbreak)
    data <- with(fake_outbreak, outbreaker_data(dates = collecDates, w_dens = w, dna = dat$dna))
    config <- outbreaker_config(data = data)
    ll <- outbreaker2:::create_loglike(data)
    param <- outbreaker_create_mcmc(data = data, config = config)$current

    ## compute likelihoods
    out <- ll$timing(param)

    ## test expected values
    expect_is(out, "numeric")
    expect_equal(out, -162.133443661423)

    ## test that likelihoods add up
    expect_equal(out, ll$timing_sampling(param) + ll$timing_infections(param))

})




## test ll$all ##
test_that("ll$all gives expected results", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    data(fake_outbreak)
    data <- with(fake_outbreak, outbreaker_data(dates = collecDates, w_dens = w, dna = dat$dna))
    config <- outbreaker_config(data = data)
    ll <- outbreaker2:::create_loglike(data)
    param <- outbreaker_create_mcmc(data = data, config = config)$current

    ## compute likelihoods
    out <- ll$all(param = param)
    out_timing <- ll$timing(param = param)
    out_genetic <- ll$genetic(param = param)
    out_reporting <- ll$reporting(param = param)

    ## test expected values
    expect_is(out, "numeric")
    expect_equal(out, -1115.21438541)

    ## test that likelihoods add up
    expect_equal(out_timing + out_genetic + out_reporting, out)

})







## test local likelihoods ##
test_that("ll$all with i specified gives expected results", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    data(fake_outbreak)
    data <- with(fake_outbreak, outbreaker_data(dates = collecDates, w_dens = w, dna = dat$dna))
    config <- outbreaker_config(data = data)
    ll <- outbreaker2:::create_loglike(data)
    param <- outbreaker_create_mcmc(data = data, config = config)$current

    ## compute local likelihoods
    sum_local_timing_sampling <- sum(sapply(seq_len(data$N), ll$timing_sampling, param = param))
    sum_local_timing_infections <- sum(sapply(seq_len(data$N), ll$timing_infections, param = param))
    sum_local_timing <- sum(sapply(seq_len(data$N), ll$timing, param = param))
    sum_local_genetic <- sum(sapply(seq_len(data$N), ll$genetic, param = param))
    sum_local_reporting <- sum(sapply(seq_len(data$N), ll$reporting, param = param))
    sum_local_all <- sum(sapply(seq_len(data$N), ll$all, param = param))

    out_timing <- ll$timing(param = param)
    out_timing_sampling <- ll$timing_sampling(param = param)
    out_timing_infections <- ll$timing_infections(param = param)
    out_genetic <- ll$genetic(param = param)
    out_reporting <- ll$reporting(param = param)
    out_all <- ll$all(param = param)

    ## tests sum of local against global
    expect_equal(sum_local_timing_sampling, out_timing_sampling)
    expect_equal(sum_local_timing_infections, out_timing_infections)
    expect_equal(sum_local_timing, out_timing)
    expect_equal(sum_local_genetic, out_genetic)
    expect_equal(sum_local_reporting, out_reporting)
    expect_equal(sum_local_all, out_all)

    ## test internal sums add up
    expect_equal(sum_local_timing_sampling + sum_local_timing_infections, sum_local_timing)
    expect_equal(sum_local_timing + sum_local_genetic + sum_local_reporting, sum_local_all)

})




## test create_loglike ##
test_that("create_loglike create functions with closure", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    data <- outbreaker_data()
    out <- outbreaker2:::create_loglike(data)

    ## tests
    expect_is(out, "list")
    expect_equal(length(out), 6)
    expect_equal(names(out), c("genetic",
                               "reporting",
                               "timing_infections",
                               "timing_sampling",
                               "timing",
                               "all"
                              )
                 )

    ## check that all items are functions
    expect_true(all(vapply(out, is.function, logical(1))))

    ## check that closure worked
    expect_identical(data, environment(out$genetic)$data)
    expect_identical(data, environment(out$reporting)$data)
    expect_identical(data, environment(out$timing_infections)$data)
    expect_identical(data, environment(out$timing_sampling)$data)

})






test_that("likelihood functions return -Inf when needed", {
    ## skip on CRAN
    skip_on_cran()

    ## generate data
    times <- 4:0
    alpha <- c(NA,rep(1,4))
    w <- c(.1, .2, .5, .2, .1)
    data <- outbreaker_data(dates = times, w_dens = w)
    config <- outbreaker_config(data = data, init_tree = alpha)
    ll <- outbreaker2:::create_loglike(data)
    param <- outbreaker_create_mcmc(data = data, config = config)$current
    few_cases <- as.integer(c(1,3,4))
    rnd_cases <- sample(sample(seq_len(data$N), 3, replace = FALSE))


    ## test ll$timing_infection ##
    out_infections <- ll$timing_infections(param)
    out_infections_few_cases <- ll$timing_infections(param, few_cases)
    out_infections_rnd_cases <- ll$timing_infections(param, rnd_cases)
    ref <- .ll_timing_infections(data, param)
    ref_few_cases <- .ll_timing_infections(data, param, few_cases)
    ref_rnd_cases <- .ll_timing_infections(data, param, rnd_cases)

    ## test values
    expect_is(out_infections, "numeric")
    expect_equal(out_infections, -Inf)
    expect_equal(out_infections_few_cases, -Inf)

    ## test against reference
    expect_equal(out_infections, ref)
    expect_equal(out_infections_few_cases, ref_few_cases)
    expect_equal(out_infections_rnd_cases, ref_rnd_cases)




    ## test ll$timing_sampling ##
    old_t_inf <- param$t_inf
    param$t_inf <- times
    out_sampling <- ll$timing_sampling(param)
    out_sampling_few_cases <- ll$timing_sampling(param, few_cases)
    out_sampling_rnd_cases <- ll$timing_sampling(param, rnd_cases)
    ref <- .ll_timing_sampling(data, param)
    ref_few_cases <- .ll_timing_sampling(data, param, few_cases)
    ref_rnd_cases <- .ll_timing_sampling(data, param, rnd_cases)
    param$t_inf <- old_t_inf

    ## test values
    expect_is(out_sampling, "numeric")
    expect_equal(out_sampling, -Inf)
    expect_equal(out_sampling_few_cases, -Inf)

    ## test against reference
    expect_equal(out_sampling, ref)
    expect_equal(out_sampling_few_cases, ref_few_cases)
    expect_equal(out_sampling_rnd_cases, ref_rnd_cases)




    ## test ll$timing ##
    out_timing <- ll$timing(param)
    out_timing_few_cases <- ll$timing(param, few_cases)
    out_timing_rnd_cases <- ll$timing(param, rnd_cases)

    ## test values
    expect_is(out_timing, "numeric")
    expect_equal(out_timing, -Inf)
    expect_equal(out_timing_few_cases, -Inf)



    ## test ll$all ##
    out_all <- ll$all(param)
    out_all_few_cases <- ll$all(param, few_cases)
    out_all_rnd_cases <- ll$all(param, rnd_cases)

    ## test values
    expect_is(out_all, "numeric")
    expect_equal(out_all, -Inf)
    expect_equal(out_all_few_cases, -Inf)

    ## test against reference
    expect_equal(out_all, ref)
    expect_equal(out_all_few_cases, ref_few_cases)
    expect_equal(out_all_rnd_cases, ref_rnd_cases)



})




