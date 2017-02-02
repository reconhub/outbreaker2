context("Test likelihood functions")




test_that("ll_timing_infections gives expected results", {
    ## skip on CRAN
    skip_on_cran()

    ## generate data
    times <- 0:4
    alpha <- c(NA,rep(1,4))
    w <- c(.1, .2, .5, .2, .1)
    data <- outbreaker_data(dates = times, w_dens = w)
    config <- create_config(data = data, init_tree = alpha)
    param <- create_param(data = data, config = config)$current
    few_cases <- as.integer(c(1,3,4))
    rnd_cases <- sample(sample(seq_len(data$N), 3, replace = FALSE))


    ## tests
    out <- cpp_ll_timing_infections(data, param)
    out_few_cases <- cpp_ll_timing_infections(data, param, few_cases)
    out_rnd_cases <- cpp_ll_timing_infections(data, param, rnd_cases)
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




## test cpp_ll_timing_sampling ##
test_that("cpp_ll_timing_sampling gives expected results", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    times <- 0:4
    alpha <- c(NA,rep(1,4))
    samp_times <- times + c(1, 1, 2, 3, 4)
    f <- c(.1, .2, .5, .2, .1)
    data <- outbreaker_data(dates = samp_times, f_dens = f)
    config <- create_config(data = data, init_t_inf = times, init_tree = alpha)
    param <- create_param(data = data, config = config)$current
    few_cases <- as.integer(c(1,3,4))
    rnd_cases <- sample(sample(seq_len(data$N), 3, replace = FALSE))


    ## tests
    out <- cpp_ll_timing_sampling(data, param)
    out_few_cases <- cpp_ll_timing_sampling(data, param, few_cases)
    out_rnd_cases <- cpp_ll_timing_sampling(data, param, rnd_cases)
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




## test cpp_ll_genetic ##
test_that("cpp_ll_genetic gives expected results", {
    ## skip on CRAN ##
    skip_on_cran()

    ## generate data ##
    data(fake_outbreak)
    data <- with(fake_outbreak,
                 outbreaker_data(dates = sample,
                                 w_dens = w,
                                 dna = dna))
    config <- create_config(data = data, init_mu = 0.543e-4)
    param <- create_param(data = data, config = config)$current
    few_cases <- as.integer(c(1,3,4))
    rnd_cases <- sample(sample(seq_len(data$N), 5, replace = FALSE))

    ## tests ##
    ## expected values
    out <- cpp_ll_genetic(data, param)
    out_few_cases <- cpp_ll_genetic(data, param, few_cases)
    out_rnd_cases <- cpp_ll_genetic(data, param, rnd_cases)
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






## test cpp_ll_reporting ##
test_that("cpp_ll_reporting gives expected results", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    data(fake_outbreak)
    data <- with(fake_outbreak,
                 outbreaker_data(dates = sample, w_dens = w, dna = dna))
    config <- create_config(data = data, init_mu = 0.543e-4)
    param <- create_param(data = data, config = config)$current
    few_cases <- as.integer(c(1,3,4))
    rnd_cases <- sample(sample(seq_len(data$N), 3, replace = FALSE))


    ## tests
    out <- cpp_ll_reporting(data, param)
    out_few_cases <- cpp_ll_reporting(data, param, few_cases)
    out_rnd_cases <- cpp_ll_reporting(data, param, rnd_cases)
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






## test cpp_ll_timing ##
test_that("cpp_ll_timing gives expected results", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    data(fake_outbreak)
    data <- with(fake_outbreak,
                 outbreaker_data(dates = sample, w_dens = w, dna = dna))
    config <- create_config(data = data)
    param <- create_param(data = data, config = config)$current

    ## compute likelihoods
    out <- cpp_ll_timing(data, param)

    ## test expected values
    expect_is(out, "numeric")
    expect_equal(out, -162.133443661423)

    ## test that likelihoods add up
    expect_equal(out, cpp_ll_timing_sampling(data, param) + cpp_ll_timing_infections(data, param))

})




## test cpp_ll_all ##
test_that("cpp_ll_all gives expected results", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    data(fake_outbreak)
    data <- with(fake_outbreak,
                 outbreaker_data(dates = sample, w_dens = w, dna = dna))
    config <- create_config(data = data)
    param <- create_param(data = data, config = config)$current

    ## compute likelihoods
    out <- cpp_ll_all(data, param = param)
    out_timing <- cpp_ll_timing(data, param = param)
    out_genetic <- cpp_ll_genetic(data, param = param)
    out_reporting <- cpp_ll_reporting(data, param = param)

    ## test expected values
    expect_is(out, "numeric")
    expect_equal(out, -1115.21438541)

    ## test that likelihoods add up
    expect_equal(out_timing + out_genetic + out_reporting, out)

})







## test local likelihoods ##
test_that("cpp_ll_all with i specified gives expected results", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    data(fake_outbreak)
    data <- with(fake_outbreak,
                 outbreaker_data(dates = sample, w_dens = w, dna = dna))
    config <- create_config(data = data)
    param <- create_param(data = data, config = config)$current

    ## compute local likelihoods
    sum_local_timing_sampling <- sum(sapply(seq_len(data$N),
                                            function(i) cpp_ll_timing_sampling(data, param, i)))
    
    sum_local_timing_infections <- sum(sapply(seq_len(data$N), 
                                            function(i) cpp_ll_timing_infections(data, param, i)))
    
    sum_local_timing <- sum(sapply(seq_len(data$N), 
                                            function(i) cpp_ll_timing(data, param, i)))
    
    sum_local_genetic <- sum(sapply(seq_len(data$N), 
                                            function(i) cpp_ll_genetic(data, param, i)))
    
    sum_local_reporting <- sum(sapply(seq_len(data$N), 
                                            function(i) cpp_ll_reporting(data, param, i)))
    
    sum_local_all <- sum(sapply(seq_len(data$N), 
                                            function(i) cpp_ll_all(data, param, i)))

    out_timing <- cpp_ll_timing(data, param = param)
    out_timing_sampling <- cpp_ll_timing_sampling(data, param = param)
    out_timing_infections <- cpp_ll_timing_infections(data, param = param)
    out_genetic <- cpp_ll_genetic(data, param = param)
    out_reporting <- cpp_ll_reporting(data, param = param)
    out_all <- cpp_ll_all(data, param = param)

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





test_that("likelihood functions return -Inf when needed", {
    ## skip on CRAN
    skip_on_cran()

    ## generate data
    times <- 4:0
    alpha <- c(NA,rep(1,4))
    w <- c(.1, .2, .5, .2, .1)
    data <- outbreaker_data(dates = times, w_dens = w)
    config <- create_config(data = data, init_tree = alpha)
    param <- create_param(data = data, config = config)$current
    few_cases <- as.integer(c(1,3,4))
    rnd_cases <- sample(sample(seq_len(data$N), 3, replace = FALSE))


    ## test cpp_ll_timing_infection ##
    out_infections <- cpp_ll_timing_infections(data, param)
    out_infections_few_cases <- cpp_ll_timing_infections(data, param, few_cases)
    out_infections_rnd_cases <- cpp_ll_timing_infections(data, param, rnd_cases)
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




    ## test cpp_ll_timing_sampling ##
    old_t_inf <- param$t_inf
    param$t_inf <- times
    out_sampling <- cpp_ll_timing_sampling(data, param)
    out_sampling_few_cases <- cpp_ll_timing_sampling(data, param, few_cases)
    out_sampling_rnd_cases <- cpp_ll_timing_sampling(data, param, rnd_cases)
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




    ## test cpp_ll_timing ##
    out_timing <- cpp_ll_timing(data, param)
    out_timing_few_cases <- cpp_ll_timing(data, param, few_cases)
    out_timing_rnd_cases <- cpp_ll_timing(data, param, rnd_cases)

    ## test values
    expect_is(out_timing, "numeric")
    expect_equal(out_timing, -Inf)
    expect_equal(out_timing_few_cases, -Inf)



    ## test cpp_ll_all ##
    out_all <- cpp_ll_all(data, param)
    out_all_few_cases <- cpp_ll_all(data, param, few_cases)
    out_all_rnd_cases <- cpp_ll_all(data, param, rnd_cases)

    ## test values
    expect_is(out_all, "numeric")
    expect_equal(out_all, -Inf)
    expect_equal(out_all_few_cases, -Inf)

    ## test against reference
    expect_equal(out_all, ref)
    expect_equal(out_all_few_cases, ref_few_cases)
    expect_equal(out_all_rnd_cases, ref_rnd_cases)

})






test_that("Customisation with identical functions works", {

    ## skip on CRAN
    skip_on_cran()

    ## check custom_likelihoods
    expect_identical(custom_likelihoods(),
                     custom_likelihoods(custom_likelihoods()))
    
    ## generate data
    data(fake_outbreak)
    data <- with(fake_outbreak,
                 outbreaker_data(dates = sample,
                                 w_dens = w,
                                 dna = dna))
    config <- create_config(data = data, init_mu = 0.543e-4)
    param <- create_param(data = data, config = config)$current
    few_cases <- as.integer(c(1,3,4))
    rnd_cases <- sample(sample(seq_len(data$N), 5, replace = FALSE))

  
    ## generate custom functions with 2 arguments
    f_genetic <- function(data, param) cpp_ll_genetic(data, param)
    f_timing_infections  <-  function(data, param) cpp_ll_timing_infections(data, param)
    f_timing_sampling  <-  function(data, param) cpp_ll_timing_sampling(data, param)
    f_reporting  <-  function(data, param) cpp_ll_reporting(data, param)
    
    list_functions <- custom_likelihoods(genetic = f_genetic,
                       timing_infections = f_timing_infections,
                       timing_sampling = f_timing_sampling,
                       reporting = f_reporting)


    ## tests
    expect_equal(cpp_ll_genetic(data, param),
                 cpp_ll_genetic(data, param, , f_genetic))

    expect_equal(cpp_ll_timing_infections(data, param),
                 cpp_ll_timing_infections(data, param, , f_timing_infections))
    

    expect_equal(cpp_ll_timing_sampling(data, param),
                 cpp_ll_timing_sampling(data, param, , f_timing_sampling))
    
    expect_equal(cpp_ll_timing_sampling(data, param),
                 cpp_ll_timing_sampling(data, param, , f_timing_sampling))
    
    expect_equal(cpp_ll_reporting(data, param),
                 cpp_ll_reporting(data, param, , f_reporting))
    
    expect_equal(cpp_ll_timing(data, param),
                 cpp_ll_timing(data, param, , list_functions))

    expect_equal(cpp_ll_all(data, param),
                 cpp_ll_all(data, param, , list_functions))
    
})












test_that("Customisation with pi-returning functions works", {

    ## skip on CRAN
    skip_on_cran()

    ## generate data ##
    data(fake_outbreak)
    data <- with(fake_outbreak,
                 outbreaker_data(dates = sample,
                                 w_dens = w,
                                 dna = dna))
    config <- create_config(data = data, init_mu = 0.543e-4)
    param <- create_param(data = data, config = config)$current
    few_cases <- as.integer(c(1,3,4))
    rnd_cases <- sample(sample(seq_len(data$N), 5, replace = FALSE))

  
    ## generate custom functions with 2 arguments
    f <- function(data, param) return(pi);
    
    list_functions <- custom_likelihoods(genetic = f,
                       timing_infections = f,
                       timing_sampling = f,
                       reporting = f)


    ## tests
    expect_equal(pi,
                 cpp_ll_genetic(data, param, , f))

    expect_equal(pi, 
                 cpp_ll_timing_infections(data, param, , f))
    
    expect_equal(pi,
                 cpp_ll_timing_sampling(data, param, , f))
    
    expect_equal(pi,
                 cpp_ll_timing_sampling(data, param, , f))
    
    expect_equal(pi,
                 cpp_ll_reporting(data, param, , f))
    
    expect_equal(2 * pi,
                 cpp_ll_timing(data, param, , list_functions))

    expect_equal(4 * pi,
                 cpp_ll_all(data, param, , list_functions))
    
})
