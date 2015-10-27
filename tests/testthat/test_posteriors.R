context("Test posterior functions")


## test ll.timing.infections ##
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
    like <- ll.all(times=times, ances=ances, log.w=log(w), D=D, mu=mu, gen.length=gen.length)
    prior <- prior.all(mu=mu)
    post <- post.all(times=times, D=D, gen.length=gen.length, log.w=log(w), ances=ances, mu=mu)
    post.check <- like+prior
    expect_equal(post, post.check)
})
