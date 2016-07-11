context("Test likelihood functions")


## test ll.timing.infections ##
test_that("ll.timing.infections gives expected results", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    times <- 0:4
    alpha <- c(NA,rep(1,4))
    w <- c(.1, .2, .5, .2, .1)
    data <- outbreaker.data(dates=times, w.dens=w)
    config <- outbreaker.config(data=data, init.tree=alpha)
    ll <- outbreaker2:::create.loglike(data)
    param <- outbreaker.create.mcmc(data=data, config=config)

    ## tests
    out <- ll$timing.infections(param)
    expect_is(out, "numeric")
    expect_equal(out, -6.214608098)
})




## test ll$timing.sampling ##
test_that("ll$timing.sampling gives expected results", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    times <- 0:4
    alpha <- c(NA,rep(1,4))
    samp.times <- times + c(1, 1, 2, 3, 4)
    f <- c(.1, .2, .5, .2, .1)
    data <- outbreaker.data(dates=samp.times, f.dens=f)
    ll <- outbreaker2:::create.loglike(data)
    config <- outbreaker.config(data=data, init.t.inf=times, init.tree=alpha)
    param <- outbreaker.create.mcmc(data=data, config=config)

    ## tests
    out <- ll$timing.sampling(param=param)
    expect_is(out, "numeric")
    expect_equal(out, -8.51719319142)
})




## test ll$genetic ##
test_that("ll$genetic gives expected results", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    data(fake.outbreak)
    data <- with(fake.outbreak, outbreaker.data(dates=collecDates, w.dens=w, dna=dat$dna))
    ll <- outbreaker2:::create.loglike(data)
    config <- outbreaker.config(data=data, init.mu=0.543e-4)
    param <- outbreaker.create.mcmc(data=data, config=config)
    fewcases <- as.integer(c(1,3,4))


    ## tests
    out <- ll$genetic(param)
    out.fewcases <- ll$genetic(param, fewcases)
    expect_is(out, "numeric")
    expect_equal(out, -997.840630502)
    expect_equal(out.fewcases, -266.251194283819)

    expect_equal(ll_genetic(data, param, integer(0)), -997.840630502)
    expect_equal(ll_genetic(data, param, fewcases), -266.251194283819)

})






## test ll$reporting ##
test_that("ll$reporting gives expected results", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    data(fake.outbreak)
    data <- with(fake.outbreak, outbreaker.data(dates=collecDates, w.dens=w, dna=dat$dna))
    config <- outbreaker.config(data=data, init.mu=0.543e-4)
    ll <- outbreaker2:::create.loglike(data)
    param <- outbreaker.create.mcmc(data=data, config=config)


    ## tests
    out <- ll$reporting(param=param)
    expect_is(out, "numeric")
    expect_equal(out, -3.05545495408)

})





## test ll$all ##
test_that("ll$all gives expected results", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    data(fake.outbreak)
    data <- with(fake.outbreak, outbreaker.data(dates=collecDates, w.dens=w, dna=dat$dna))
    config <- outbreaker.config(data=data)
    ll <- outbreaker2:::create.loglike(data)
    param <- outbreaker.create.mcmc(data=data, config=config)

    ## compute likelihoods
    out <- ll$all(param=param)
    out.timing <- ll$timing(param=param)
    out.genetic <- ll$genetic(param=param)
    out.reporting <- ll$reporting(param=param)

    ## test expected values
    expect_is(out, "numeric")
    expect_equal(out, -1115.21438541)

    ## test that likelihoods add up
    expect_equal(out.timing + out.genetic + out.reporting, out)

})







## test local likelihoods ##
test_that("ll$all with i specified gives expected results", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    data(fake.outbreak)
    data <- with(fake.outbreak, outbreaker.data(dates=collecDates, w.dens=w, dna=dat$dna))
    config <- outbreaker.config(data=data)
    ll <- outbreaker2:::create.loglike(data)
    param <- outbreaker.create.mcmc(data=data, config=config)

    ## compute local likelihoods
    sum.local.timing.sampling <- sum(sapply(seq_len(data$N), ll$timing.sampling, param=param))
    sum.local.timing.infections <- sum(sapply(seq_len(data$N), ll$timing.infections, param=param))
    sum.local.timing <- sum(sapply(seq_len(data$N), ll$timing, param=param))
    sum.local.genetic <- sum(sapply(seq_len(data$N), ll$genetic, param=param))
    sum.local.reporting <- sum(sapply(seq_len(data$N), ll$reporting, param=param))
    sum.local.all <- sum(sapply(seq_len(data$N), ll$all, param=param))

    out.timing <- ll$timing(param=param)
    out.timing.sampling <- ll$timing.sampling(param=param)
    out.timing.infections <- ll$timing.infections(param=param)
    out.genetic <- ll$genetic(param=param)
    out.reporting <- ll$reporting(param=param)
    out.all <- ll$all(param=param)

    ## tests sum of local against global
    expect_equal(sum.local.timing.sampling, out.timing.sampling)
    expect_equal(sum.local.timing.infections, out.timing.infections)
    expect_equal(sum.local.timing, out.timing)
    expect_equal(sum.local.genetic, out.genetic)
    expect_equal(sum.local.reporting, out.reporting)
    expect_equal(sum.local.all, out.all)

    ## test internal sums add up
    expect_equal(sum.local.timing.sampling + sum.local.timing.infections, sum.local.timing)
    expect_equal(sum.local.timing + sum.local.genetic + sum.local.reporting, sum.local.all)

})




## test create.loglike ##
test_that("create.loglike create functions with closure", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    data <- outbreaker.data()
    out <- outbreaker2:::create.loglike(data)

    ## tests
    expect_is(out, "list")
    expect_equal(length(out), 6)
    expect_equal(names(out), c("genetic",
                               "reporting",
                               "timing.infections",
                               "timing.sampling",
                               "timing",
                               "all"
                              )
                 )

    ## check that all items are functions
    expect_true(all(vapply(out, is.function, logical(1))))

    ## check that closure worked
    expect_identical(data, environment(out$genetic)$data)
    expect_identical(data, environment(out$reporting)$data)
    expect_identical(data, environment(out$timing.infections)$data)
    expect_identical(data, environment(out$timing.sampling)$data)

})
