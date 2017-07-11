context("Test outbreaker")

## test output format ##
test_that("outbreaker's output have expected format", {
    ## skip on CRAN
    skip_on_cran()


    ## get data
    data(fake_outbreak)
    x <- fake_outbreak

    ## run outbreaker
    data <- list(dna = x$dna, dates = x$onset, w_dens = x$w)
    config <- list(n_iter = 10, sample_every = 1, paranoid = TRUE,
                   find_import = FALSE)
    out <- outbreaker(data, config)


    out_df <- as.data.frame(out)

    data <- list(dates = x$onset, w_dens = x$w)
    out2 <- outbreaker(data, config)

    ## check output
    expect_is(out, "outbreaker_chains")
    expect_is(out_df, "data.frame")
    expect_equal(nrow(out), 10)
    expect_true(!any(is.na(out_df$post)))
    expect_true(all(out_df$post> -1e30))

})




## test convergence in various settings ##
test_that("results ok: DNA + time", {
    ## skip on CRAN
    skip_on_cran()


    ## get data
    data(fake_outbreak)
    x <- fake_outbreak

    ## outbreaker DNA + time ##
    ## analysis
    set.seed(1)
    data <- list(dna = x$dna, dates = x$onset, w_dens = x$w)
    config <- list(n_iter = 5000, sample_every = 50,
                   init_tree = "seqTrack", find_import = TRUE,
                   move_pi = FALSE)

    out <- outbreaker(data = data, config = config)

    ## checks
    out_smry <- summary(out, burnin = 1000)

    ## approx log post values
    expect_true(min(out_smry$post) > -1250)

    ## at least 80% ancestries correct
    temp <- mean(out_smry$tree$from==x$ances, na.rm = TRUE)
    expect_true(temp > .8)

    ## infection date within 3 days on average
    temp <- mean(abs(out_smry$tree$time - x$onset), na.rm = TRUE)
    expect_true(temp < 3.5)

    ## mu within reasonable ranges
    expect_true(out_smry$mu[2] > 1e-5 &&
                out_smry$mu[5] < 2e-4)

})





test_that("results ok: time, no DNA", {
    ## skip on CRAN
    skip_on_cran()


    ## get data
    data(fake_outbreak)
    x <- fake_outbreak

    ## outbreaker time, no DNA ##
    ## analysis
    set.seed(1)

    data <- list(dates = x$onset, w_dens = x$w)
    config <- list(n_iter = 5000, sample_every = 50,
                   init_tree="star", find_import = FALSE,
                   move_kappa = FALSE)

    out_no_dna <- outbreaker(data = data, config = config)

    ## checks
    out_no_dna.smry <- summary(out_no_dna, burnin = 1000)

    ## approx log post values
    expect_true(min(out_no_dna.smry$post) > -90)

    ## infection date within 3 days on average
    temp <- mean(abs(out_no_dna.smry$tree$time - x$onset), na.rm = TRUE)
    expect_true(temp < 3.5)

    ## check that support for ancestries is weak
    sup <- na.omit(out_no_dna.smry$tree$support)
    expect_lt(quantile(sup, .9), .5)
    expect_lt(mean(sup), .35)

})



test_that("results ok: time, ctd, DNA", {
    ## skip on CRAN
    skip_on_cran()


    ## get data
    data(fake_outbreak)
    x <- fake_outbreak

    ## outbreaker time, ctd, no DNA ##
    ## analysis
    set.seed(1)

    tTree <- data.frame(i = x$ances,
                        j = 1:length(x$ances))
    
    ctd <- sim_ctd(tTree, eps = 0.9, lambda = 0.1)

    data <- list(dates = x$onset, w_dens = x$w,
                 ctd = ctd[,1:2], dna = x$dna)
    config <- list(n_iter = 5000, sample_every = 50,
                   init_tree="star", find_import = FALSE,
                   move_kappa = FALSE)

    out_ctd <- outbreaker(data = data, config = config)

    ## checks
    out_ctd.smry <- summary(out_ctd, burnin = 1000)

    ## at least 90% ancestries correct
    temp <- mean(out_ctd.smry$tree$from==x$ances, na.rm = TRUE)
    expect_true(temp > .9)

    ## eps and lambda parameter estimates within reasonable ranges
    quant_eps <- quantile(out_ctd$eps, c(0.25, 0.75))
    expect_true(quant_eps[[1]] > 0.7 &
                quant_eps[[2]] < 0.95)

    quant_lambda <- quantile(out_ctd$lambda, c(0.25, 0.75))
    expect_true(quant_lambda[[1]] > 0.05 &
                quant_lambda[[2]] < 0.35)
        
    ## approx log post values
    expect_true(min(out_ctd.smry$post) > -1200)

})




test_that("results ok: easy run, no missing cases, no import", {
    ## skip on CRAN
    skip_on_cran()


    ## get data
    data(fake_outbreak)
    x <- fake_outbreak

    ## outbreaker, no missing cases ##
    ## analysis
    set.seed(1)
    config <-  list(n_iter = 5e3, sample_every = 50,
                    init_tree = "seqTrack", move_kappa = FALSE,
                    move_pi = FALSE, init_pi = 1,
                    find_import = FALSE)
    data <- list(dna = x$dna, dates = x$onset, w_dens = x$w)

    out_no_missing <- outbreaker(data = data, config = config)

    ## checks
    out_no_missing_smry <- summary(out_no_missing, burnin = 3500)

    ## approx log post values
    expect_true(min(out_no_missing_smry$post) > -935)

    ## at least 85% ancestries correct
    temp <- mean(out_no_missing_smry$tree$from==x$ances, na.rm = TRUE)
    expect_true(temp > .85)

    ## infection datewithin 3 days on average
    temp <- mean(abs(out_no_missing_smry$tree$time - x$onset), na.rm = TRUE)
    expect_true(temp < 3.5)

})



test_that("results ok: no missing cases, detect imported cases",
{
    ## skip on CRAN
    skip_on_cran()


    ## get data
    data(fake_outbreak)
    x <- fake_outbreak

    ## outbreaker, no missing cases, detect imported cases ##
    ## analysis
    set.seed(1)
    out_with_import <- outbreaker(data = list(dna = x$dna, dates = x$onset,
                                              w_dens = x$w),
                                  config = list(n_iter = 5000, sample_every = 100,
                                                init_tree="star",
                                                move_kappa = FALSE, move_pi = FALSE,
                                                init_pi = 1, find_import = TRUE)
                                  )

    ## checks
    out_with_import_smry <- summary(out_with_import, burnin = 500)
    out_with_import_smry$tree$from[is.na(out_with_import_smry$tree$from)] <- 0
    x$ances[is.na(x$ances)] <- 0

    ## approx log post values
    expect_true(min(out_with_import_smry$post) > -460)

    ## at least 80% ancestries correct
    temp <- mean(out_with_import_smry$tree$from==x$ances, na.rm = TRUE)
    expect_true(temp >= .80)

    ## infection datewithin 3 days on average
    temp <- mean(abs(out_with_import_smry$tree$time - x$onset), na.rm = TRUE)
    expect_true(temp < 3.5)

})






test_that("results ok: kappa and pi", {
    ## skip on CRAN
    skip_on_cran()


    ## get data
    onset <- c(0, 2, 6, 13)
    w <- c(.25, .5, .25)

    ## outbreaker, no missing cases, detect imported cases ##
    ## analysis
    set.seed(1)

    data <- list(dates = onset, w_dens = w)
    config <- list(n_iter = 5000, sample_every = 50, init_tree = "star",
                   move_kappa = TRUE, move_pi = TRUE, init_pi = 1,
                   find_import = FALSE, max_kappa = 10)


    out <- outbreaker(data, config)

    smry <- summary(out, burnin = 500)

    ## checks
    expect_equal(smry$tree$from, c(NA, 1, 2, 3))
    expect_equal(smry$tree$generations, c(NA, 1, 2, 4.5))
    expect_true(min(smry$post) > -35)
    expect_true(all(smry$pi[3:4] > 0.5 & smry$pi[3:4] < 0.8))

})







test_that("results ok: outbreaker with custom priors",
{
    ## skip on CRAN
    skip_on_cran()


    ## get data
    data(fake_outbreak)
    x <- fake_outbreak

    data <- list(dna = x$dna, dates = x$onset, w_dens = x$w)

    config <-  list(n_iter = 5e3, sample_every = 50,
                    init_tree = "star", move_kappa = TRUE,
                    move_pi = TRUE, init_pi = 1,
                    find_import = FALSE)


    ## plot(function(x) dbeta(x, 100, 1)) # to see the distribution

    f_pi_1 <- function(x) dbeta(x$pi, 100, 1, log = TRUE) # stronger prior
    f_pi_2 <- function(x) 0.0 # flat prior

    set.seed(1)

    out_1 <- outbreaker(data, config, priors = list(pi = f_pi_1))
    out_2 <- outbreaker(data, config, priors = list(pi = f_pi_2))

    smry_1 <- summary(out_1, burnin = 500)
    smry_2 <- summary(out_2, burnin = 500)

    expect_true(smry_1$pi[2] > 0.9)
    expect_true(smry_2$pi[2] > 0.2 && smry_2$pi[5] < 0.6)

})








test_that("results ok: outbreaker with fixed number returning priors and likelihoods",
{
    ## skip on CRAN
    skip_on_cran()


    ## get data

    data(fake_outbreak)
    x <- fake_outbreak

    data <- list(dna = x$dna, dates = x$onset, w_dens = x$w)

    config <-  list(n_iter = 1000, sample_every = 10,
                    init_tree = "star", move_kappa = TRUE,
                    move_pi = TRUE, init_pi = 1,
                    find_import = FALSE)


    ## custom functions

    f1 <- function(x) return(0.0)
    f2 <- function(x, y) return(0.0)
    f3 <- function(x) return(1.123)

    priors1 <- custom_priors(mu = f1, pi = f1)
    priors2 <- custom_priors(mu = f1, pi = f3)

    ll <- custom_likelihoods(genetic = f2,
                             timing_infections = f2,
                             timing_sampling = f2,
                             reporting = f2)


    ## tests
    out1 <- outbreaker(data, config, priors1, ll)
    out2 <- outbreaker(data, config, priors2, ll)

    expect_true(all(out1$post == 0))
    expect_true(all(out2$post == 1.123))

})
