context("Test non-exported functions")


## test .choose_possible_alpha ##
test_that("test: .choose_possible_alpha", {
  ## skip on CRAN
  skip_on_cran()


  ## get data
  ans1 <- outbreaker2:::.choose_possible_alpha(1:10, 1)
  ans2 <- outbreaker2:::.choose_possible_alpha(1:10, 10)
  ans3 <- outbreaker2:::.choose_possible_alpha(1:10, 2)
  ans4 <- outbreaker2:::.are_possible_alpha(1:10, 5)
  ans5 <- outbreaker2:::.are_possible_alpha(c(1:4, 1:3, 1), 2)
  ans6 <- outbreaker2:::.are_possible_alpha(1:10, 1)

  ## check output
  expect_true(is.na(ans1))
  expect_is(ans2, "integer")
  expect_true(ans2<10 || ans2>0)
  expect_equal(ans3, 1L)
  expect_error(.choose_possible_alpha(1:10, NA),
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
  t_inf <- c(2, 1, 3, 4, 4, 3)
  data <- outbreaker_data(dates = t_inf+1)
  config <- create_config(init_tree = alpha, init_t_inf = t_inf, data = data)
  config2 <- config
  config2$move_alpha <- c(rep(TRUE,4),FALSE, TRUE)
  param <- create_param(config = config, data = data)$current

  ## test can be alpha
  expect_equal(can_move_alpha(param, config), c(TRUE, FALSE,rep(TRUE,4)))
  expect_equal(can_move_alpha(param, config2), c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE))

  ## test ancestor selection
  set.seed(1)
  to_move <- replicate(10, select_alpha_to_move(param,config))
  expect_equal(to_move, c(3,3,4,6,3,6,6,5,5,1))
  to_move2 <- replicate(10, select_alpha_to_move(param,config2))
  expect_equal(to_move2, c(1,1,4,3,6,3,4,6,3,6))

})






## Test swapping of ancestries in Cpp
test_that("Ancestries swapping in Rcpp works", {
  skip_on_cran()

  ## make dummy tree
  param_old <- list(alpha = c(NA, 1L, 2L, 3L, 2L),
                    t_inf = c(1L, 2L, 3L, 4L, 3L),
                    kappa = c(1L, 2L, 1L, 1L, 1L))

  ## no swapping exception: i is imported
  param_new <- cpp_swap_cases(param_old, 1L)
  expect_identical(param_old, param_new)

  ## swap even when ancestor is imported
  param_new <- cpp_swap_cases(param_old, 2L)
  expect_equal(param_new$alpha, c(2L, NA, 1L, 3L, 1L))
  expect_equal(param_new$t_inf, c(2L, 1L, 3L, 4L, 3L))
  expect_equal(param_new$kappa, c(2L, 1L, 1L, 1L, 1L))

  ## swap 3 and 2
  param_new <- cpp_swap_cases(param_old, 3L)
  expect_equal(param_new$alpha, c(NA, 3L, 1L, 2L, 3L))
  expect_equal(param_new$t_inf, c(1L, 3L, 2L, 4L, 3L))

})






## test find_descendents ##
test_that("Testing find_descendents", {
  ## skip on CRAN
  skip_on_cran()


  ## generate data
  data <- with(fake_outbreak,
               outbreaker_data(dates = sample, w_dens = w, dna = dna))

  config <- create_config(data = data)
  param <- create_param(data = data, config = config)$current

  ## tests
  expect_equal(find_descendents(param, 1), c(2,4,28))
  expect_equal(find_descendents(param, 30), integer(0))

  ## test cpp version
  expect_equal(cpp_find_descendents(c(NA), 10), integer(0))
  expect_equal(cpp_find_descendents(c(1,NA,2,1), 1), c(1L, 4L))
  expect_equal(cpp_find_descendents(c(NA, 1,1,1,2,2,NA,1,1), 1),
               c(2L, 3L, 4L, 8L, 9L))
  expect_equal(cpp_find_descendents(c(NA,1,NA,2), NA_integer_), integer(0))
})






test_that("Testing add_convolutions", {
  ## skip on CRAN
  skip_on_cran()


  ## generate data
  data(fake_outbreak)
  data <- with(fake_outbreak,
               outbreaker_data(dates = sample,
                               w_dens = w,
                               dna = dna))
  config <- create_config(data = data)

  out <- add_convolutions(data = data, config = config)

  ## tests
  expect_is(out$log_w_dens, "matrix")
  expect_equal(dim(out$log_w_dens), c(5, length(out$w_dens)))
  expect_true(!any(is.na(out$log_w_dens)))

})






test_that("Testing cpp_find_local_cases", {
  ## skip on CRAN
  skip_on_cran()

  expect_equal(cpp_find_local_cases(NA, 1), 1L)
  expect_equal(cpp_find_local_cases(c(NA, 1, 2), 1), c(1L, 2L))
  expect_equal(cpp_find_local_cases(c(NA, 1, 1, 1), 1),
               as.integer(1:4))
  expect_equal(sort(cpp_find_local_cases(c(NA, 1, 1, 1), 2)),
               as.integer(1:4))

  tre <- c(NA, 1, 1, 3, 3, 5, 2)
  expect_equal(cpp_find_local_cases(tre, 1), 1:3)
  expect_equal(sort(cpp_find_local_cases(tre, 2)),
               c(1, 2, 3, 7))
  expect_equal(sort(cpp_find_local_cases(tre, 7)),
               c(2, 7))
  expect_equal(sort(cpp_find_local_cases(tre, 4)),
               c(3, 4, 5))
  expect_equal(sort(cpp_find_local_cases(tre, 5)),
               c(3, 4, 5, 6))
  expect_equal(sort(cpp_find_local_cases(tre, 6)),
               c(5, 6))

})






test_that("Testing cpp_get_n_mutations", {
  ## skip on CRAN
  skip_on_cran()


  ## sequences in natural order

  data <- outbreaker_data(dates = fake_outbreak$onset,
                          dna = fake_outbreak$dna)
  D <- data$D

  expect_equal(cpp_get_n_mutations(data, 1, 1), 0L)
  expect_equal(cpp_get_n_mutations(data, 11, 3), D[11, 3])
  expect_equal(cpp_get_n_mutations(data, 21, 30), D[21, 30])


  ## sequences in reverse order, named

  dna.rev <- fake_outbreak$dna
  rownames(dna.rev) <- as.character(30:1)

  data.rev <- outbreaker_data(dates = fake_outbreak$onset,
                              dna = dna.rev)
  D.rev <- data.rev$D

  expect_equal(rownames(D.rev), rownames(dna.rev))
  expect_equal(cpp_get_n_mutations(data.rev, 11, 3), D.rev["11", "3"])
  expect_equal(cpp_get_n_mutations(data.rev, 21, 30), D.rev["21", "30"])


  ## with some sequences
  dna.miss <- fake_outbreak$dna[c(3, 4, 11), ]
  rownames(dna.miss) <- c("3", "4", "11")

  data.miss <- outbreaker_data(dates = fake_outbreak$onset,
                               dna = dna.miss)
  D.miss <- data.miss$D

  expect_equal(rownames(D.miss), rownames(dna.miss))
  expect_equal(cpp_get_n_mutations(data.miss, 3, 4), D.miss["3", "4"])
  expect_equal(cpp_get_n_mutations(data.miss, 3, 4), D[3, 4])
  expect_equal(cpp_get_n_mutations(data.miss, 11, 3), D.miss["11", "3"])

  msg <- "Trying to get genetic distances between missing sequences."
  expect_error(cpp_get_n_mutations(data.miss, 1, 2), msg)


})






test_that("Look for sequenced ancestors", {

  ## skip on CRAN
  skip_on_cran()


  ## make inputs
  onset <- as.integer(c(0, 1, 2, 3, 2, 1, 2, 3))
  data <- outbreaker_data(onset = onset)

  has_dna <- as.logical(c(0, 1, 0, 1, 1, 0, 0, 1))
  data$has_dna <- has_dna

  alpha <- as.integer(c(NA, 1, 2, 3, 2, 1, 6, 7))
  kappa <- as.integer(c(NA, 1, 1, 2, 1, 2, 1, 1))
  param <- list(alpha = alpha,
                kappa = kappa)


  ## make tests
  res <- lapply(as.integer(1:8),
                function(i) cpp_lookup_sequenced_ancestor(data, param, i))

  res <- Reduce(rbind, lapply(res, data.frame))
  exp_alpha <- c(NA, NA, NA, 2L, 2L, NA, NA, NA)
  exp_n_generations <- c(NA, NA, NA, 3L, 1L, NA, NA, NA)
  exp_found <- c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE)

  expect_equal(res$alpha, exp_alpha)
  expect_equal(res$n_generations, exp_n_generations)
  expect_equal(res$found_sequenced_ancestor, exp_found)

})


