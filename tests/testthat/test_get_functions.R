context("Test get_* functions")

## shared setup
get_short_chain <- function() {
  data(fake_outbreak)
  x <- fake_outbreak
  data <- list(dna = x$dna, dates = x$onset, w_dens = x$w)
  config <- list(n_iter = 50, sample_every = 5, find_import = FALSE)
  outbreaker(data, config)
}


test_that("get_trees", {
  skip_on_cran()
  out <- get_short_chain()

  trees <- get_trees(out)
  expect_is(trees, "list")
  expect_equal(length(trees), nrow(out))
  expect_true(all(c("from", "to") %in% names(trees[[1]])))
  expect_equal(nrow(trees[[1]]), 30L)

  ## optional columns
  data(fake_outbreak)
  trees2 <- get_trees(
    out,
    kappa = TRUE,
    t_inf = TRUE,
    date = fake_outbreak$onset
  )
  expect_true(all(
    c(
      "from_kappa",
      "to_kappa",
      "from_t_inf",
      "to_t_inf",
      "from_date",
      "to_date"
    ) %in%
      names(trees2[[1]])
  ))
})


test_that("get_Ri", {
  skip_on_cran()
  out <- get_short_chain()

  ri <- get_Ri(out)
  expect_true(all(c("case", "mean", "lwr", "upr") %in% names(ri)))
  expect_equal(nrow(ri), 30L)
  expect_true(all(ri$mean >= 0))

  ri_raw <- get_Ri(out, raw = TRUE)
  expect_true(all(c("step", "case", "ri") %in% names(ri_raw)))
  expect_equal(nrow(ri_raw), nrow(out) * 30L)
})


test_that("get_offspring", {
  skip_on_cran()
  out <- get_short_chain()

  off <- get_offspring(out)
  expect_true(all(c("x", "mean", "lwr", "upr") %in% names(off)))
  expect_equal(min(off$x), 0L)
  expect_equal(sum(off$mean), 1, tolerance = 1e-10)

  off_raw <- get_offspring(out, raw = TRUE)
  expect_true(all(c("step", "x", "y") %in% names(off_raw)))
  pmf_sums <- tapply(off_raw$y, off_raw$step, sum)
  expect_true(all(abs(pmf_sums - 1) < 1e-10))
})


test_that("get_entropy", {
  skip_on_cran()
  out <- get_short_chain()

  ent <- get_entropy(out)
  expect_equal(length(ent), 30L)
  expect_true(all(ent >= 0 & ent <= 1))

  ent_raw <- get_entropy(out, normalise = FALSE)
  expect_true(all(ent_raw >= 0))
})


test_that("get_si", {
  skip_on_cran()
  data(fake_outbreak)
  x <- fake_outbreak
  dat <- outbreaker_data(dates = x$onset, dna = x$dna, w_dens = x$w)
  config <- list(n_iter = 50, sample_every = 5, find_import = FALSE)
  out <- outbreaker(dat, config)

  si <- get_si(out, dat)
  expect_true(all(c("x", "mean", "lwr", "upr") %in% names(si)))
  expect_equal(sum(si$mean), 1, tolerance = 1e-10)

  si_raw <- get_si(out, dat, raw = TRUE)
  expect_true(all(c("step", "x", "y") %in% names(si_raw)))
  pmf_sums <- tapply(si_raw$y, si_raw$step, sum)
  expect_true(all(abs(pmf_sums - 1) < 1e-10))
})


test_that("get_accuracy", {
  skip_on_cran()
  out <- get_short_chain()
  data(fake_outbreak)

  true_tree <- data.frame(
    from = as.character(fake_outbreak$ances),
    to = as.character(seq_along(fake_outbreak$onset)),
    stringsAsFactors = FALSE
  )

  acc <- get_accuracy(out, true_tree)
  expect_equal(length(acc), nrow(out))
  expect_true(all(acc >= 0 & acc <= 1))
})
