context("Test sanitize_pmf")


## test settings ##
test_that("test: settings are processed fine", {
  ## skip on CRAN
  skip_on_cran()
  
  
  ## get data
  x <- c(0,0,0.1,0.2,0.3,0.2,0,0.2,0, 0)
  
  ## check output
  expect_is(x, "numeric")
  
  expect_is(tails_pmf(x), "numeric")
  
  expect_equal(length(tails_pmf(x)), length(x))
  
  expect_equal(sum(tails_pmf(x)), 1)
  
  expect_false(0 %in% tails_pmf(x)[1])
  expect_false(0 %in% tails_pmf(x)[length(tails_pmf(x))])
  
  
  expect_error(tails_pmf(c(0,0,0,0,0)))
  
  expect_error(tails_pmf(c("x", "y", "z", 1)), "non-numeric argument to mathematical function")
  
})




