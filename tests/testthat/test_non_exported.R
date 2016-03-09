context("Test non-exported functions")


## test find.possible.ances ##
test_that("test: find.possible.ances", {
    ## skip on CRAN
    skip_on_cran()
    rm(list=ls())

    ## get data
    ans1 <- outbreaker2:::find.possible.ances(1:10, 1)
    ans2 <- outbreaker2:::find.possible.ances(1:10, 10)
    ans3 <- outbreaker2:::find.possible.ances(1:10, 2)

    ## check output
    expect_true(is.na(ans1))
    expect_is(ans2, "integer")
    expect_true(ans2<10 || ans2>0)
    expect_equal(ans3, 1L)
    expect_error(find.possible.ances(1:10, NA))
})


## test swap.ances ##
test_that("Auxiliary functions for ancestries are working", {
    ## skip on CRAN
    skip_on_cran()
    rm(list=ls())

    ## generate data
    ances <- c(NA,1,1,2,2)
    t.inf <- c(0, 2,3, 5,7)
    data <- outbreaker.data(dates=t.inf+1)
    config <- outbreaker.config(init.tree=ances, init.t.inf=t.inf)
    config2 <- config
    config2$move.ances <- c(rep(TRUE,4),FALSE)
    param <- outbreaker.mcmc.init(config=config, data=data)

    ## test can be ances
    expect_equal(can.move.ances(param, config), c(FALSE,rep(TRUE,4)))
    expect_equal(can.move.ances(param, config2), c(FALSE, TRUE, TRUE, TRUE, FALSE))

    ## test ancestor selection
    set.seed(1)
    to.move <- replicate(10,select.ances.to.move(param,config))
    expect_equal(to.move, c(3,3,4,5,2,5,5,4,4,2))
    to.move2 <- replicate(10,select.ances.to.move(param,config2))
    expect_equal(to.move2, c(2,2,4,3,4,3,4,4,3,4))

    ## tests non swapping (imported case)
    expect_warning(res <- outbreaker2:::swap.ances(param, 1))
    expect_equal(param, res)

    ## test full swapping
    res <- outbreaker2:::swap.ances(param, 2)
    expect_equal(res$current.ances, c(2,NA,2,1,1))
    expect_equal(res$current.t.inf, c(2,0,3,5,7))

})
