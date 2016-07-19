context("Test likelihood functions")



######################################################
######################################################
## End of the 'reference' functions. Tests are below.
## test ll.timing.infections ##
######################################################
######################################################

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
    param <- outbreaker.create.mcmc(data=data, config=config)$current
    few.cases <- as.integer(c(1,3,4))
    rnd.cases <- sample(sample(seq_len(data$N), 3, replace=FALSE))


    ## tests
    out <- ll$timing.infections(param)
    out.few.cases <- ll$timing.infections(param, few.cases)
    out.rnd.cases <- ll$timing.infections(param, rnd.cases)
    ref <- .ll.timing.infections(data, param)
    ref.few.cases <- .ll.timing.infections(data, param, few.cases)
    ref.rnd.cases <- .ll.timing.infections(data, param, rnd.cases)

    expect_is(out, "numeric")
    expect_equal(out, -6.214608098)
    expect_equal(out.few.cases, -2.30258509299405)

    ## test against reference
    expect_equal(out, ref)
    expect_equal(out.few.cases, ref.few.cases)
    expect_equal(out.rnd.cases, ref.rnd.cases)
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
    param <- outbreaker.create.mcmc(data=data, config=config)$current
    few.cases <- as.integer(c(1,3,4))
    rnd.cases <- sample(sample(seq_len(data$N), 3, replace=FALSE))


    ## tests
    out <- ll$timing.sampling(param)
    out.few.cases <- ll$timing.sampling(param, few.cases)
    out.rnd.cases <- ll$timing.sampling(param, rnd.cases)
    ref <- .ll.timing.sampling(data, param)
    ref.few.cases <- .ll.timing.sampling(data, param, few.cases)
    ref.rnd.cases <- .ll.timing.sampling(data, param, rnd.cases)

    expect_is(out, "numeric")
    expect_equal(out, -8.51719319142)
    expect_equal(out.few.cases, -4.60517018598809)

    ## test against reference
    expect_equal(out, ref)
    expect_equal(out.few.cases, ref.few.cases)
    expect_equal(out.rnd.cases, ref.rnd.cases)

})




## test ll$genetic ##
test_that("ll$genetic gives expected results", {
    ## skip on CRAN ##
    skip_on_cran()

    ## generate data ##
    data(fake.outbreak)
    data <- with(fake.outbreak, outbreaker.data(dates=collecDates, w.dens=w, dna=dat$dna))
    ll <- outbreaker2:::create.loglike(data)
    config <- outbreaker.config(data=data, init.mu=0.543e-4)
    param <- outbreaker.create.mcmc(data=data, config=config)$current
    few.cases <- as.integer(c(1,3,4))
    rnd.cases <- sample(sample(seq_len(data$N), 5, replace=FALSE))

    ## tests ##
    ## expected values
    out <- ll$genetic(param)
    out.few.cases <- ll$genetic(param, few.cases)
    out.rnd.cases <- ll$genetic(param, rnd.cases)
    ref <- .ll.genetic(data, param)
    ref.few.cases <- .ll.genetic(data, param, few.cases)
    ref.rnd.cases <- .ll.genetic(data, param, rnd.cases)

    expect_is(out, "numeric")
    expect_equal(out, -997.840630502)
    expect_equal(out.few.cases, -266.251194283819)

    ## test against reference
    expect_equal(out, ref)
    expect_equal(out.few.cases, ref.few.cases)
    expect_equal(out.rnd.cases, ref.rnd.cases)

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
    param <- outbreaker.create.mcmc(data=data, config=config)$current
    few.cases <- as.integer(c(1,3,4))
    rnd.cases <- sample(sample(seq_len(data$N), 3, replace=FALSE))


    ## tests
    out <- ll$reporting(param)
    out.few.cases <- ll$reporting(param, few.cases)
    out.rnd.cases <- ll$reporting(param, rnd.cases)
    ref <- .ll.reporting(data, param)
    ref.few.cases <- .ll.reporting(data, param, few.cases)
    ref.rnd.cases <- .ll.reporting(data, param, rnd.cases)

    expect_is(out, "numeric")
    expect_equal(out, -3.05545495407696)
    expect_equal(out.few.cases, -0.210721031315653)

    ## test against reference
    expect_equal(out, ref)
    expect_equal(out.few.cases, ref.few.cases)
    expect_equal(out.rnd.cases, ref.rnd.cases)

})






## test ll$timing ##
test_that("ll$timing gives expected results", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    data(fake.outbreak)
    data <- with(fake.outbreak, outbreaker.data(dates=collecDates, w.dens=w, dna=dat$dna))
    config <- outbreaker.config(data=data)
    ll <- outbreaker2:::create.loglike(data)
    param <- outbreaker.create.mcmc(data=data, config=config)$current

    ## compute likelihoods
    out <- ll$timing(param)

    ## test expected values
    expect_is(out, "numeric")
    expect_equal(out, -162.133443661423)

    ## test that likelihoods add up
    expect_equal(out, ll$timing.sampling(param) + ll$timing.infections(param))

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
    param <- outbreaker.create.mcmc(data=data, config=config)$current

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
    param <- outbreaker.create.mcmc(data=data, config=config)$current

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






test_that("likelihood functions return -Inf when needed", {
    ## skip on CRAN
    skip_on_cran()

    ## generate data
    times <- 4:0
    alpha <- c(NA,rep(1,4))
    w <- c(.1, .2, .5, .2, .1)
    data <- outbreaker.data(dates=times, w.dens=w)
    config <- outbreaker.config(data=data, init.tree=alpha)
    ll <- outbreaker2:::create.loglike(data)
    param <- outbreaker.create.mcmc(data=data, config=config)$current
    few.cases <- as.integer(c(1,3,4))
    rnd.cases <- sample(sample(seq_len(data$N), 3, replace=FALSE))


    ## test ll$timing.infection ##
    out.infections <- ll$timing.infections(param)
    out.infections.few.cases <- ll$timing.infections(param, few.cases)
    out.infections.rnd.cases <- ll$timing.infections(param, rnd.cases)
    ref <- .ll.timing.infections(data, param)
    ref.few.cases <- .ll.timing.infections(data, param, few.cases)
    ref.rnd.cases <- .ll.timing.infections(data, param, rnd.cases)

    ## test values
    expect_is(out.infections, "numeric")
    expect_equal(out.infections, -Inf)
    expect_equal(out.infections.few.cases, -Inf)

    ## test against reference
    expect_equal(out.infections, ref)
    expect_equal(out.infections.few.cases, ref.few.cases)
    expect_equal(out.infections.rnd.cases, ref.rnd.cases)




    ## test ll$timing.sampling ##
    old.t.inf <- param$t.inf
    param$t.inf <- times
    out.sampling <- ll$timing.sampling(param)
    out.sampling.few.cases <- ll$timing.sampling(param, few.cases)
    out.sampling.rnd.cases <- ll$timing.sampling(param, rnd.cases)
    ref <- .ll.timing.sampling(data, param)
    ref.few.cases <- .ll.timing.sampling(data, param, few.cases)
    ref.rnd.cases <- .ll.timing.sampling(data, param, rnd.cases)
    param$t.inf <- old.t.inf

    ## test values
    expect_is(out.sampling, "numeric")
    expect_equal(out.sampling, -Inf)
    expect_equal(out.sampling.few.cases, -Inf)

    ## test against reference
    expect_equal(out.sampling, ref)
    expect_equal(out.sampling.few.cases, ref.few.cases)
    expect_equal(out.sampling.rnd.cases, ref.rnd.cases)




    ## test ll$timing ##
    out.timing <- ll$timing(param)
    out.timing.few.cases <- ll$timing(param, few.cases)
    out.timing.rnd.cases <- ll$timing(param, rnd.cases)

    ## test values
    expect_is(out.timing, "numeric")
    expect_equal(out.timing, -Inf)
    expect_equal(out.timing.few.cases, -Inf)



    ## test ll$all ##
    out.all <- ll$all(param)
    out.all.few.cases <- ll$all(param, few.cases)
    out.all.rnd.cases <- ll$all(param, rnd.cases)

    ## test values
    expect_is(out.all, "numeric")
    expect_equal(out.all, -Inf)
    expect_equal(out.all.few.cases, -Inf)

    ## test against reference
    expect_equal(out.all, ref)
    expect_equal(out.all.few.cases, ref.few.cases)
    expect_equal(out.all.rnd.cases, ref.rnd.cases)



})




