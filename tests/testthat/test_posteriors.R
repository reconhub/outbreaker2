context("Test posterior functions")


## test posterior computation ##
test_that("post = like + prior", {
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
    like <- ll.all(t.inf=times, sampling.times=times+2, ances=ances, log.w=log(w), log.f=log(w), D=D, mu=mu, gen.length=gen.length)
    prior <- prior.all(mu=mu)
    post <- post.all(t.inf=times, sampling.times=times+2, D=D, gen.length=gen.length, log.w=log(w), log.f=log(w), ances=ances, mu=mu)
    post.check <- like+prior
    expect_equal(post, post.check)
})



## test posterior computation ##
test_that("posterior is atomic", {
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
    like <- ll.all(t.inf=times, sampling.times=times+2, ances=ances, log.w=log(w), log.f=log(w), D=D, mu=mu, gen.length=gen.length)
    prior <- prior.all(mu=mu)
    post <- post.all(t.inf=times, sampling.times=times+2, D=D, gen.length=gen.length, log.w=log(w), log.f=log(w), ances=ances, mu=mu)

    expect_equal(length(like), 1)
    expect_equal(length(prior), 1)
    expect_equal(length(post), 1)
})
