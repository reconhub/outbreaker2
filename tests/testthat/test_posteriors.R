context("Test posterior functions")


## test create.posteriors ##
test_that("create.posteriors gives expected results", {
    ## skip on CRAN
    skip_on_cran()


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
    expect_identical(ll$genetic, environment(out$genetic)$loglike$genetic)
    expect_identical(ll$reporting, environment(out$reporting)$loglike$reporting)
    expect_identical(priors$mu, environment(out$genetic)$prior.mu)
    expect_identical(priors$pi, environment(out$reporting)$prior.pi)

})





## test posterior computation ##
test_that("posteriors computations give expected values", {
    ## skip on CRAN
    skip_on_cran()


    ## generate inputs
    data(fake.outbreak)
    data <- with(fake.outbreak, outbreaker.data(dates=collecDates, w.dens=w, dna=dat$dna))
    config <- outbreaker.config(data=data)
    param <- outbreaker.create.mcmc(data=data, config=config)
    ll <- create.loglike(data)
    priors <- create.priors(config)
    posteriors <- create.posteriors(ll, priors)

    ## tests
    expect_equal(posteriors$genetic(param), ll$genetic(param) + priors$mu(param))
    expect_equal(posteriors$reporting(param), ll$reporting(param) + priors$pi(param))
    expect_equal(posteriors$all(param), ll$all(param) + priors$all(param))
})

