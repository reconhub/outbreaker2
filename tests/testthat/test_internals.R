context("Test internal functions")


## test ll.timing.infections ##
test_that("ll.timing.infections gives expected results", {
    ## skip on CRAN
    skip_on_cran()

    ## generate data
    times <- 0:4
    ances <- c(NA,rep(1,4))
    w <- c(0, .1, .2, .5, .2, .1)

    ## tests
    out <- ll.timing.infections(times=times, ances=ances, log.w=log(w))
    expect_is(out, "numeric")
    expect_equal(out, -6.214608098)
})




## test ll.timing.sampling ##
test_that("ll.timing.sampling gives expected results", {
    ## skip on CRAN
    skip_on_cran()

    ## generate data
    times <- 0:4
    samp.times <- times + c(1, 1, 2, 3, 4)
    f <- c(0, .1, .2, .5, .2, .1)

    ## tests
    out <- ll.timing.sampling(times=times, sampling.times=samp.times, log.f=log(f))
    expect_is(out, "numeric")
    expect_equal(out, -8.51719319142)
})




## test ll.genetic ##
test_that("ll.genetic gives expected results", {
    ## skip on CRAN
    skip_on_cran()

    ## generate data
    set.seed(1)
    D <- as.matrix(round(dist(rnorm(5))))
    ances <- c(NA, 1, 1, 1, 2)
    mu <- 0.543e-4
    gen.length <- 2e4

    ## tests
    out <- ll.genetic(D=D, ances=ances, mu=mu, gen.length=gen.length)
    expect_is(out, "numeric")
    expect_equal(out, -33.80691403)
})




## test ll.all ##
test_that("ll.all gives expected results", {
    ## skip on CRAN
    skip_on_cran()

    ## generate data
    set.seed(1)
    times <- 0:4
    w <- c(0, .1, .2, .5, .2, .1)
    D <- as.matrix(round(dist(rnorm(5))))
    ances <- c(NA, 1, 1, 1, 2)
    mu <- 0.543e-4
    gen.length <- 2e4

    ## tests
    out <- ll.all(times=times, ances=ances, log.w=log(w), D=D, mu=mu, gen.length=gen.length)
    expect_is(out, "numeric")
    expect_equal(out, -39.1052313991)
})


