context("Test likelihood functions")


## test ll.timing.infections ##
test_that("ll.timing.infections gives expected results", {
    ## skip on CRAN
    skip_on_cran()
    rm(list=ls())

    ## generate data
    times <- 0:4
    ances <- c(NA,rep(1,4))
    w <- c(.1, .2, .5, .2, .1)
    data <- outbreaker.data(dates=times, w.dens=w)
    config <- outbreaker.config(data=data, init.tree=ances)
    param <- outbreaker.mcmc.init(data=data, config=config)

    ## tests
    out <- ll.timing.infections(data=data, param=param)
    expect_is(out, "numeric")
    expect_equal(out, -6.214608098)
})




## test ll.timing.sampling ##
test_that("ll.timing.sampling gives expected results", {
    ## skip on CRAN
    skip_on_cran()
    rm(list=ls())

    ## generate data
    times <- 0:4
    ances <- c(NA,rep(1,4))
    samp.times <- times + c(1, 1, 2, 3, 4)
    f <- c(.1, .2, .5, .2, .1)
    data <- outbreaker.data(dates=samp.times, f.dens=f)
    config <- outbreaker.config(data=data, init.t.inf=times, init.tree=ances)
    param <- outbreaker.mcmc.init(data=data, config=config)

    ## tests
    out <- ll.timing.sampling(data=data, param=param)
    expect_is(out, "numeric")
    expect_equal(out, -8.51719319142)
})




## test ll.genetic ##
test_that("ll.genetic gives expected results", {
    ## skip on CRAN
    skip_on_cran()
    rm(list=ls())

    ## generate data
    data(fakeOutbreak)
    data <- with(fakeOutbreak, outbreaker.data(dates=collecDates, w.dens=w, dna=dat$dna))
    config <- outbreaker.config(data=data, init.mu=0.543e-4)
    param <- outbreaker.mcmc.init(data=data, config=config)


    ## tests
    out <- ll.genetic(data=data, param=param)
    expect_is(out, "numeric")
    expect_equal(out, -997.840630502)
})




## test ll.all ##
test_that("ll.all gives expected results", {
    ## skip on CRAN
    skip_on_cran()
    rm(list=ls())

    ## generate data
    data(fakeOutbreak)
    data <- with(fakeOutbreak, outbreaker.data(dates=collecDates, w.dens=w, dna=dat$dna))
    config <- outbreaker.config(data=data)
    param <- outbreaker.mcmc.init(data=data, config=config)

    ## compute likelihoods
    out <- ll.all(data=data, param=param)
    out.timing <- ll.timing(data=data, param=param)
    out.genetic <- ll.genetic(data=data, param=param)

    ## test expected values
    expect_is(out, "numeric")
    expect_equal(out, -1112.15893046)

    ## test that likelihoods add up
    expect_equal(out.timing + out.genetic, out)
})







## test local likelihoods ##
test_that("ll.all.i gives expected results", {
    ## skip on CRAN
    skip_on_cran()
    rm(list=ls())

    ## generate data
    data(fakeOutbreak)
    data <- with(fakeOutbreak, outbreaker.data(dates=collecDates, w.dens=w, dna=dat$dna))
    config <- outbreaker.config(data=data)
    param <- outbreaker.mcmc.init(data=data, config=config)

    ## compute local likelihoods
    sum.all.timing.sampling <- sum(sapply(seq.int(data$N), ll.timing.sampling.i, data=data, param=param))
    sum.all.timing.infections <- sum(sapply(seq.int(data$N), ll.timing.infections.i, data=data, param=param))
    sum.all.timing <- sum(sapply(seq.int(data$N), ll.timing.i, data=data, param=param))
    sum.all.genetic <- sum(sapply(seq.int(data$N), ll.genetic.i, data=data, param=param))
    sum.all <- sum(sapply(seq.int(data$N), ll.all.i, data=data, param=param))

    out.timing <- ll.timing(data=data, param=param)
    out.timing.sampling <- ll.timing.sampling(data=data, param=param)
    out.timing.infections <- ll.timing.infections(data=data, param=param)
    out.genetic <- ll.genetic(data=data, param=param)
    out.all <- ll.all(data=data, param=param)

    ## tests sum of local against global
    expect_equal(sum.all.timing.sampling, out.timing.sampling)
    expect_equal(sum.all.timing.infections, out.timing.infections)
    expect_equal(sum.all.timing, out.timing)
    expect_equal(sum.all.genetic, out.genetic)
    expect_equal(sum.all, out.all)

    ## test internal sums add up
    expect_equal(sum.all.timing.sampling + sum.timing.infections, sum.all.timing)
    expect_equal(sum.all.timing + sum.all.genetic, sum.all)

})




## test outbreaker.create.loglike ##
test_that("outbreaker.create.loglike gives expected results", {
    ## skip on CRAN
    skip_on_cran()
    rm(list=ls())

    ## generate data
    out <- outbreaker.create.loglike()

    ## tests
    expect_is(out, "list")
    expect_equal(length(out), 10)
    expect_equal(names(out), c("timing.infections",
                               "timing.sampling",
                               "timing",
                               "genetic",
                               "all",
                               "timing.infections.i",
                               "timing.sampling.i",
                               "timing.i",
                               "genetic.i",
                               "all.i")
                 )

})
