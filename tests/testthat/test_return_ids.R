context("Test return_ids")

test_that("return_ids = TRUE produces same alpha values as FALSE", {
  skip_on_cran()

  data(fake_outbreak)

  o2_data <- outbreaker_data(
    dates = as.Date(fake_outbreak$sample, origin = "2020-01-01"),
    ids = paste0("Case_", 1:30),
    w_dens = fake_outbreak$w
  )

  ## run with return_ids = FALSE (default)
  set.seed(123)
  o2_noid <- outbreaker(
    data = o2_data,
    config = list(n_iter = 500, sample_every = 50, return_ids = FALSE)
  )

  ## run with return_ids = TRUE
  set.seed(123)
  o2_id <- outbreaker(
    data = o2_data,
    config = list(n_iter = 500, sample_every = 50, return_ids = TRUE)
  )

  ## extract alpha columns
  alpha_noid <- o2_noid[, grep("alpha", names(o2_noid))]
  alpha_id <- o2_id[, grep("alpha", names(o2_id))]

  ## strip "Case_" from values to get back to numeric strings
  alpha_id_stripped <- as.data.frame(
    lapply(alpha_id, function(col) {
      x <- as.integer(sub("Case_", "", col))
      x[is.na(col)] <- NA_integer_
      x
    })
  )

  ## strip "Case_" from column names
  names(alpha_id_stripped) <- sub("Case_", "", names(alpha_id_stripped))

  ## column names should match (alpha_1, alpha_2, ...)
  expect_identical(names(alpha_noid), names(alpha_id_stripped))

  ## values should be identical
  expect_identical(
    as.matrix(alpha_noid),
    as.matrix(alpha_id_stripped)
  )
})
