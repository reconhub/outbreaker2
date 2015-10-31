context("Test posterior functions")


## test posterior computation ##
test_that("post = like + prior", {
    ## skip on CRAN
    skip_on_cran()

    ## generate data
    data(fakeOutbreak)
    dat <- with(fakeOutbreak, outbreaker.data(dates=collecDates, w.dens=w, dna=dat$dna))
    config <- outbreaker.config(data=dat)
    chain <- outbreaker.mcmc.init(data=dat, config=config)

    ## tests
    like <-ll.all(data=dat, chain=chain)
    prior <- prior.all(chain)
    post <- post.all(data=dat, chain=chain)
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
    chain <- outbreaker.mcmc.init(data=dat, config=config)

    ## tests
    like <-ll.all(data=dat, chain=chain)
    prior <- prior.all(chain)
    post <- post.all(data=dat, chain=chain)

    expect_equal(length(like), 1)
    expect_equal(length(prior), 1)
    expect_equal(length(post), 1)
})
