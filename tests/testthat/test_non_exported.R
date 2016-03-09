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
    ances <- c(2, NA, 1, 3, 3, 1)
    t.inf <- c(2, 1, 3, 4, 4, 3)
    data <- outbreaker.data(dates=t.inf+1)
    config <- outbreaker.config(init.tree=ances, init.t.inf=t.inf, data=data)
    config2 <- config
    config2$move.ances <- c(rep(TRUE,4),FALSE, TRUE)
    param <- outbreaker.mcmc.init(config=config, data=data)

    ## test can be ances
    expect_equal(can.move.ances(param, config), c(TRUE, FALSE,rep(TRUE,4)))
    expect_equal(can.move.ances(param, config2), c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE))

    ## test ancestor selection
    set.seed(1)
    to.move <- replicate(10,select.ances.to.move(param,config))
    expect_equal(to.move, c(3,3,4,6,3,6,6,5,5,1))
    to.move2 <- replicate(10,select.ances.to.move(param,config2))
    expect_equal(to.move2, c(1,1,4,3,6,3,4,6,3,6))

    ## tests non swapping (imported case)
    expect_warning(res <- outbreaker2:::swap.ances(param, config, 2))
    expect_equal(param, res)

    ## trying swap x->i to i->x when x is imported (forbidden)
    res <- outbreaker2:::swap.ances(param, config, 1)
    expect_equal(res$current.ances, param$current.ances)

    ## test full swapping
    res <- outbreaker2:::swap.ances(param, config, 3)
    expect_equal(res$current.ances, c(3,NA,2,1,1,3))
    expect_equal(res$current.t.inf, c(3,1,2,4,4,3))

})
