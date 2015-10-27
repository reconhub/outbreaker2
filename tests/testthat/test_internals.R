context("Test internal functions")

times <- 0:4
ances <- c(NA,rep(1,4))
w <- c(0, .1, .2, .5, .2, .1 )

test_that("ll.timing.infections gives expected results", {
  skip_on_cran()
  out <- ll.timing.infections(times=times, ances=ances, log.w=log(w))
  expect_is(out, "numeric")
  expect_equal(out, -6.214608098)
})
