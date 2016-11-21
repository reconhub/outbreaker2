context("Test non-exported functions")


## test .choose.possible.alpha ##
test_that("test: .choose.possible.alpha", {
    ## skip on CRAN
    skip_on_cran()


    ## get data
    ans1 <- outbreaker2:::.choose.possible.alpha(1:10, 1)
    ans2 <- outbreaker2:::.choose.possible.alpha(1:10, 10)
    ans3 <- outbreaker2:::.choose.possible.alpha(1:10, 2)
    ans4 <- outbreaker2:::.are.possible.alpha(1:10, 5)
    ans5 <- outbreaker2:::.are.possible.alpha(c(1:4, 1:3, 1), 2)
    ans6 <- outbreaker2:::.are.possible.alpha(1:10, 1)

    ## check output
    expect_true(is.na(ans1))
    expect_is(ans2, "integer")
    expect_true(ans2<10 || ans2>0)
    expect_equal(ans3, 1L)
    expect_error(.choose.possible.alpha(1:10, NA),
                 "missing value where TRUE/FALSE needed")
    expect_equal(ans4, 1:4)
    expect_equal(ans5, c(1,5,8))
    expect_equal(ans6, NA)
})


## ancestry-related functions (R functions)
test_that("Auxiliary functions for ancestries are working", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    alpha <- c(2, NA, 1, 3, 3, 1)
    t.inf <- c(2, 1, 3, 4, 4, 3)
    data <- outbreaker.data(dates=t.inf+1)
    config <- outbreaker.config(init.tree=alpha, init.t.inf=t.inf, data=data)
    config2 <- config
    config2$move.alpha <- c(rep(TRUE,4),FALSE, TRUE)
    param <- outbreaker.create.mcmc(config=config, data=data)$current

    ## test can be alpha
    expect_equal(can.move.alpha(param, config), c(TRUE, FALSE,rep(TRUE,4)))
    expect_equal(can.move.alpha(param, config2), c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE))

    ## test ancestor selection
    set.seed(1)
    to.move <- replicate(10,select.alpha.to.move(param,config))
    expect_equal(to.move, c(3,3,4,6,3,6,6,5,5,1))
    to.move2 <- replicate(10,select.alpha.to.move(param,config2))
    expect_equal(to.move2, c(1,1,4,3,6,3,4,6,3,6))

    ## tests non swapping (imported case)
    expect_warning(res <- outbreaker2:::swap.cases(param, config, 2))
    expect_equal(param, res)

    ## trying swap x->i to i->x when x is imported (forbidden)
    res <- outbreaker2:::swap.cases(param, config, 1)
    expect_equal(res$alpha, param$alpha)

    ## test full swapping
    res <- outbreaker2:::swap.cases(param, config, 3)
    expect_equal(res$alpha, c(3,NA,2,1,1,3))
    expect_equal(res$t.inf, c(3,1,2,4,4,3))

})




## Test swapping of ancestries in Cpp
test_that("Ancestries swapping in Rcpp works", {
    skip_on_cran()

    ## make dummy tree
    param.old <- list(alpha = c(NA, 1L, 2L, 3L, 2L),
                      t.inf = c(1L, 2L, 3L, 4L, 3L))
    
    ## no swapping exception: i is imported
    param.new <- cpp.swap.cases(param.old, 1L)
    expect_identical(param.old, param.new)

    ## no swapping exception: i's ancestor is imported
    param.new <- cpp.swap.cases(param.old, 2L)
    expect_identical(param.old, param.new)

    ## swap 3 and 2
    param.new <- cpp.swap.cases(param.old, 3L)
    expect_equal(param.new$alpha, c(NA, 3L, 1L, 2L, 3L))
    expect_equal(param.new$t.inf, c(1L, 3L, 2L, 4L, 3L))
    
})

## ## test check.i ##
## test_that("Testing check.i", {
##     ## skip on CRAN
##     skip_on_cran()
##

##     ## generate data
##     alpha <- c(2, NA, 1, 3, 3, 1)
##     t.inf <- c(2, 1, 3, 4, 4, 3)
##     data <- outbreaker.data(dates=t.inf+1)

##     ## tests
##     expect_error(check.i(data, 0:10))
##     expect_error(check.i(data, -13))
##     expect_error(check.i(data, 5:11))
##     expect_error(check.i(data, c(NA,2)))
##     expect_error(check.i(data, "1"))
##     expect_error(check.i(data, TRUE))
##     expect_equal(check.i(data, NULL), 1:6)
##     expect_equal(check.i(data, 6),6)
##     expect_equal(check.i(data, c(1,3,2)), c(1,3,2))
##     expect_equal(check.i(data, 1:6), 1:6)

## })




## test find.descendents ##
test_that("Testing find.descendents", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    data <- with(fake.outbreak,
                 outbreaker.data(dates=collecDates, w.dens=w, dna=dat$dna))
    
    config <- outbreaker.config(data=data)
    param <- outbreaker.create.mcmc(data=data, config=config)$current

    ## tests
    expect_equal(find.descendents(param, 1), c(2,4,28))
    expect_equal(find.descendents(param, 30), integer(0))

    ## test cpp version
    expect_equal(cpp.find.descendents(c(NA), 10), integer(0))
    expect_equal(cpp.find.descendents(c(1,NA,2,1), 1), c(1L, 4L))
    expect_equal(cpp.find.descendents(c(NA, 1,1,1,2,2,NA,1,1), 1),
                                      c(2L, 3L, 4L, 8L, 9L))
    expect_equal(cpp.find.descendents(c(NA,1,NA,2), NA_integer_), integer(0))
})




test_that("Testing add.convolutions", {
    ## skip on CRAN
    skip_on_cran()


    ## generate data
    data(fake.outbreak)
    data <- with(fake.outbreak, outbreaker.data(dates=collecDates, w.dens=w, dna=dat$dna))
    config <- outbreaker.config(data=data)

    out <- add.convolutions(data=data, config=config)

    ## tests
    expect_is(out$log.w.dens, "matrix")
    expect_equal(dim(out$log.w.dens), c(5, length(out$w.dens)))
    expect_true(!any(is.na(out$log.w.dens)))

})




test_that("Testing cpp.find.local.cases", {
    ## skip on CRAN
    skip_on_cran()

    expect_equal(cpp.find.local.cases(NA, 1), 1L)
    expect_equal(cpp.find.local.cases(c(NA, 1, 2), 1), c(1L, 2L))
    expect_equal(cpp.find.local.cases(c(NA, 1, 1, 1), 1),
                 as.integer(1:4))
    expect_equal(sort(cpp.find.local.cases(c(NA, 1, 1, 1), 2)), 
                 as.integer(1:4))

    tre <- c(NA, 1, 1, 3, 3, 5, 2)
    expect_equal(cpp.find.local.cases(tre, 1), 1:3)
    expect_equal(sort(cpp.find.local.cases(tre, 2)),
                 c(1, 2, 3, 7))
    expect_equal(sort(cpp.find.local.cases(tre, 7)),
                 c(2, 7))
    expect_equal(sort(cpp.find.local.cases(tre, 4)),
                 c(3, 4, 5))
    expect_equal(sort(cpp.find.local.cases(tre, 5)),
                 c(3, 4, 5, 6))
    expect_equal(sort(cpp.find.local.cases(tre, 6)),
                 c(5, 6))
      
})
