context("Test make.fast functions")


## test plausible values ##
test_that("Outputs have expected size and types", {
    ## skip on CRAN
    skip_on_cran()

    ## generate data ##
    size <- as.integer(sample(10:100,1))
    x <- make.fast.log.runif(batch.size=size)
    x1 <- make.fast.log.runif1(batch.size=size)

    ## tests ##
    ## test length of pre-generated vector
    expect_true(length(environment(x)$values)==size)
    expect_true(length(environment(x1)$values)==size)

    ## test length of ouput
    expect_true(length(x(69))==69L)
    expect_true(length(x1())==1L)

    ## test values
    set.seed(1);a <- make.fast.rnorm(mean=1, sd=1.123)(10)
    set.seed(1);b <- rnorm(10, mean=1, sd=1.123)
    expect_equal(a,b)

    ## test counters
    y1 <- make.fast.rnorm1()
    i <- environment(y1)$counter
    expect_equal(i, 0L)
    y1()
    i <- environment(y1)$counter
    expect_equal(i, 1L)
    y <- make.fast.rnorm()
    y(4)
    expect_equal(environment(y)$counter, 4L)

})
