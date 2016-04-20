context("Test posterior functions")




## test create.posteriors ##
test_that("create.posteriors gives expected results", {
    ## skip on CRAN
    skip_on_cran()
    rm(list=ls())

    ## generate data
    data <- outbreaker.data()
    config <- outbreaker.config()
    param <- outbreaker.create.mcmc(data=data, config=config)
    ll <- create.loglike(data)
    priors <- create.priors(config)
    out <- create.posteriors(ll, priors)

    ## tests
    expect_is(out, "list")
    expect_equal(length(out), 3)
    expect_equal(names(out), c("genetic",
                               "reporting",
                               "all")
                 )

    ## check that all items are functions
    expect_true(all(vapply(out, is.function, logical(1))))

    ## check that closure worked
    expect_identical(ll$genetic, environment(out$genetic)$ll.genetic)
    expect_identical(ll$reporting, environment(out$reporting)$ll.reporting)
    expect_identical(priors$mu, environment(out$genetic)$prior.mu)
    expect_identical(priors$pi, environment(out$reporting)$prior.pi)

})



## test posterior computation ##
test_that("post = like + prior", {
    ## skip on CRAN
    skip_on_cran()
    rm(list=ls())

    ## generate inputs
    data(fakeOutbreak)
    data <- with(fakeOutbreak, outbreaker.data(dates=collecDates, w.dens=w, dna=dat$dna))
    config <- outbreaker.config(data=data)
    param <- outbreaker.create.mcmc(data=data, config=config)
    ll <- create.loglike(data)


    ## generate outputs

    ## tests
    like <-ll.all(data=dat, param=param)
    prior <- prior.all(param)
    post <- post.all(data=dat, param=param)
    post.check <- like+prior
    expect_equal(post, post.check)
})



## test posterior computation ##
test_that("posterior is atomic", {
   ## skip on CRAN
    skip_on_cran()

    ## generate data
    data(fakeOutbreak)
    dat <- with(fakeOutbreak, outbreaker.data(dates=collecDates, w.dens=w, dna=dat$dna))
    config <- outbreaker.config(data=dat)
    param <- outbreaker.mcmc.init(data=dat, config=config)

    ## tests
    like <-ll.all(data=dat, param=param)
    prior <- prior.all(param)
    post <- post.all(data=dat, param=param)

    expect_equal(length(like), 1)
    expect_equal(length(prior), 1)
    expect_equal(length(post), 1)
})


