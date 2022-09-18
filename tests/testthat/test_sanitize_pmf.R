context("Test sanitize_pmf")


## test settings ##
test_that("test: settings are processed fine", {
  ## skip on CRAN
  skip_on_cran()
  
  
  ## get data
  x <- c(0,0,0.1,0.2,0.3,0.2,0,0.2,0, 0)
  
  ## check output
  expect_is(x, "numeric")
  
  expect_is(sanitize_pmf(x), "numeric")
  
  expect_equal(length(sanitize_pmf(x)), length(x))
  
  expect_equal(sum(sanitize_pmf(x)), 1)
  
  expect_false(0 %in% sanitize_pmf(x))

  expect_warning(sanitize_pmf(c(0,0,0,0,0)))
  
  expect_error(sanitize_pmf(c("x", "y", "z", 1)), "non-numeric argument to mathematical function")
  
})




