context("Test binding functions")

test_that("Function binding throws errors as expected", {
    ## skip on CRAN
    skip_on_cran()


    ## get data
    data(fake_outbreak)
    x <- fake_outbreak

    ## get data
    data <- list(dna = x$dna, dates = x$onset, w_dens = x$w)
    config <- list(n_iter = 10, sample_every = 1, paranoid = TRUE,
                   find_import = FALSE)


    ## checks errors
    msg <- "'...' is empty but 'f' has more than one argument."
    expect_error(bind_to_function(cpp_ll_genetic), msg)

    msg <- "All arguments provided through '...' need to be named."
    expect_error(bind_to_function(cpp_ll_genetic, data), msg)

    msg <- paste("Arguments of cpp_move_mu missing from '...'",
                 "with no default: config")
    expect_error(bind_to_function(cpp_move_mu, data = data),
                 msg)

})





test_that("Function binding encloses arguments", {
    ## skip on CRAN
    skip_on_cran()


    ## get data
    data(fake_outbreak)
    x <- fake_outbreak$dat

    ## get data
    data <- list(dna = x$dna, dates = x$onset, w_dens = x$w)
    config <- list(n_iter = 10, sample_every = 1, paranoid = TRUE,
                   find_import = FALSE)


    ## checks
    f <- bind_to_function(cpp_ll_genetic, data = data)
    expect_identical(environment(f)$data, data)

    f <- bind_to_function(cpp_ll_genetic, data = data, i = 1:3)
    expect_identical(environment(f)$i, 1:3)
    
})
