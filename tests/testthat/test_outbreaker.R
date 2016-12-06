context("Test outbreaker")

## test output format ##
test_that("outbreaker's output have expected format", {
    ## skip on CRAN
    skip_on_cran()


    ## get data
    data(fake.outbreak)
    dat <- fake.outbreak$dat
    w <- fake.outbreak$w

    ## run outbreaker
    out <- outbreaker(data=list(dna=dat$dna, dates=dat$onset, w.dens=w),
                      config=list(n.iter=100, sample.every=10, paranoid=TRUE))


    out.df <- as.data.frame(out)

    out2 <- outbreaker(data=list(dates=dat$onset, w.dens=w),
                       config=list(n.iter=100, sample.every=10, paranoid=TRUE))

    ## check output
    expect_is(out, "outbreaker.chains")
    expect_is(out.df, "data.frame")
    expect_equal(nrow(out), 11)
    expect_true(!any(is.na(out.df$post)))
    expect_true(all(out.df$post> -1e30))

})




## test convergence in various settings ##
test_that("decent results for toy example, DNA + time", {
    ## skip on CRAN
    skip_on_cran()


    ## get data
    data(fake.outbreak)
    dat <- fake.outbreak$dat
    w <- fake.outbreak$w

    ## outbreaker DNA + time ##
    ## analysis
    set.seed(1)
    out <- outbreaker(data=list(dna=dat$dna, dates=dat$onset, w.dens=w),
                      config=list(n.iter=5000, sample.every=100,
                      init.tree="star", find.import=FALSE))

    ## checks
    out.smry <- summary(out, burnin=1000)

    ## approx log post values
    expect_true(min(out.smry$post) > -950)

    ## at least 75% ancestries correct
    expect_true(mean(out.smry$tree$from==dat$ances, na.rm=TRUE) > .75)

    ## infection datewithin 3 days on average
    expect_true(mean(abs(out.smry$tree$time - dat$onset), na.rm=TRUE)<3.5)

    ## mu between 2e-4 and 4 e-4
    expect_true(min(out.smry$mu) > 0.0002 && max(out.smry$mu) < 0.00042)
})





test_that("decent results for toy example, time, no DNA", {
    ## skip on CRAN
    skip_on_cran()


    ## get data
    data(fake.outbreak)
    dat <- fake.outbreak$dat
    w <- fake.outbreak$w

   ## outbreaker time, no DNA ##
    ## analysis
    set.seed(1)
    out.no.dna <- outbreaker(data=list(dates=dat$onset, w.dens=w),
                             config=list(n.iter=5000, sample.every=100,
                                         init.tree="star", find.import=FALSE))

    ## checks
    out.no.dna.smry <- summary(out.no.dna, burnin=1000)
    
    ## approx log post values
    expect_true(min(out.no.dna.smry$post) > -100) 

    ## at least 5% ancestries correct
    temp <- mean(out.no.dna.smry$tree$from==dat$ances, na.rm=TRUE)
    expect_true(temp > .05)

    ## infection datewithin 3 days on average
    temp <- mean(abs(out.no.dna.smry$tree$time - dat$onset), na.rm=TRUE)
    expect_true(temp < 3.5)
})




test_that("decent results for toy example, no missing cases", {
    ## skip on CRAN
    skip_on_cran()


    ## get data
    data(fake.outbreak)
    dat <- fake.outbreak$dat
    w <- fake.outbreak$w

    ## outbreaker, no missing cases ##
    ## analysis
    set.seed(1)
    out.no.missing <- outbreaker(data=list(dna=dat$dna, dates=dat$onset, w.dens=w),
                                 config=list(n.iter=10000, sample.every=100,
                                             init.tree="star", move.kappa=FALSE,
                                             move.pi=FALSE, init.pi=1,
                                             find.import=FALSE)
                                 )

    ## checks
    out.no.missing.smry <- summary(out.no.missing, burnin=2000)

    ## approx log post values
    expect_true(min(out.no.missing.smry$post) > -935)

    ## at least 80% ancestries correct
    temp <- mean(out.no.missing.smry$tree$from==dat$ances, na.rm=TRUE)
    expect_true(temp > .8)

    ## infection datewithin 3 days on average
    temp <- mean(abs(out.no.missing.smry$tree$time - dat$onset), na.rm=TRUE)
    expect_true(temp < 3.5)
})



test_that("decent results for toy example, no missing cases, detect imported cases",
{
    ## skip on CRAN
    skip_on_cran()
    
    
    ## get data
    data(fake.outbreak)
    dat <- fake.outbreak$dat
    w <- fake.outbreak$w
    
    ## outbreaker, no missing cases, detect imported cases ##
    ## analysis
    set.seed(1)
    out.with.import <- outbreaker(data = list(dna=dat$dna, dates=dat$onset,
                                              w.dens=w),
                                  config = list(n.iter=10000, sample.every=100,
                                                init.tree="star",
                                                move.kappa=FALSE, move.pi=FALSE,
                                                init.pi=1, find.import=TRUE)
                                  )
    
    ## checks
    out.with.import.smry <- summary(out.with.import, burnin=1000)
    out.with.import.smry$tree$from[is.na(out.with.import.smry$tree$from)] <- 0
    dat$ances[is.na(dat$ances)] <- 0

    ## approx log post values
    expect_true(min(out.with.import.smry$post) > -460)

    ## at least 80% ancestries correct
    temp <- mean(out.with.import.smry$tree$from==dat$ances, na.rm=TRUE)
    expect_true(temp >= .80)

    ## infection datewithin 3 days on average
    temp <- mean(abs(out.with.import.smry$tree$time - dat$onset), na.rm=TRUE)
    expect_true(temp < 3.5) 
    
})






test_that("decent results for kappa and pi", {
    ## skip on CRAN
    skip_on_cran()


    ## get data
    onset <- c(0,2,8,20)
    w <- c(.25, .5, .25)
    
    ## outbreaker, no missing cases, detect imported cases ##
    ## analysis
    set.seed(1)

    config <- list(n.iter = 5000, sample.every = 50, init.tree = "star",
                   move.kappa = TRUE, move.pi = TRUE, init.pi = 1,
                   find.import = FALSE)
    
    out <- outbreaker(data = list(dates = onset, w.dens= w),
                                  config = config)
    ## plot(out)
    smry <- summary(out, burnin=1000)

    ## checks
    expect_equal(smry$tree$from, c(NA, 1, 2, 3))
    expect_equal(smry$tree$generations, c(NA, 1, 4, 5))
    expect_true(min(smry$post) > -28)
    expect_true(all(smry$pi[3:4] > 0.5 & smry$pi[3:4] < 0.7))

})
