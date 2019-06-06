context("Test outbreaker data and settings")


## test data ##
test_that("test: data are processed fine", {
    ## skip on CRAN
    skip_on_cran()


    ## get data
    x <- fake_outbreak
    out <- outbreaker_data(dates = x$onset, dna = x$dna, w_dens = x$w)
    out_nodna <- outbreaker_data(dates = x$onset, w_dens = x$w)

    ## check output
    expect_is(out, "list")
    expect_is(out$D, "matrix")
    expect_equal(out$max_range, 11)
    expect_equal(out_nodna$L, 0)
    expect_equal(out$L, 1e4)
    expect_equal(out$w_dens, out$f_dens)
    expect_equal(out$log_w_dens[1,], out$log_f_dens)
    expect_error(outbreaker_data(dates = 1, w_dens = c(0,-1)),
                 "w_dens has negative entries")

    expect_error(outbreaker_data(dates = 1, w_dens = c(0,1), f_dens = c(0,-1)),
                 "f_dens has negative entries")

    wrong_lab_dna <- x$dna
    rownames(wrong_lab_dna) <- paste0("host_", seq_len(nrow(wrong_lab_dna)))
    expect_error(outbreaker_data(dates = x$onset, dna = wrong_lab_dna, w_dens = x$w),
                 "DNA sequence labels don't match case ids")


})







test_that("outbreaker_data accepts epicontacts and case labelling", {
  
    ## skip on CRAN
    skip_on_cran()


    ## get data
    data(fake_outbreak)
    x <- fake_outbreak

    ## outbreaker time, ctd, no DNA ##
    ## analysis
    set.seed(1)

    ## use id labelling throughout
    ids <- replicate(length(fake_outbreak$sample),
                     paste(sample(letters, 5, TRUE), collapse = ""))

    ## make epi_contacts object
    tTree <- data.frame(i = ids[x$ances],
                        j = ids[1:length(x$ances)])
    ctd <- sim_ctd(tTree, eps = 0.9, lambda = 0.1)
    epi_c <- epicontacts::make_epicontacts(linelist = data.frame(id = ids),
                                           contacts = ctd,
                                           directed = TRUE)

    data <- outbreaker_data(dates = x$onset,
                            dna = x$dna,
                            ctd = epi_c,
                            w_dens = x$w)

    ## test recursiveness
    data <- outbreaker_data(data = data)

    ## check ids are carried through
    expect_equal(data$ids, epi_c$linelist$id)
    
    ## make sure directionality is carried through
    expect_true(data$ctd_directed)


    ## case labelling via dates
    dates <- x$onset
    names(dates) <- ids
    
    data <- outbreaker_data(dates = dates,
                            dna = x$dna,
                            ctd = ctd,
                            w_dens = x$w,
                            ctd_directed = TRUE)

    ## test recursiveness
    data <- outbreaker_data(data = data)

    ## check ids are carried through
    expect_equal(data$ids, ids)
    
    ## make sure directionality is carried through
    expect_true(data$ctd_directed)

    ## check the number of contacts are correct
    expect_equal(nrow(ctd), sum(data$contacts))


    ## toggle directionality
    data <- outbreaker_data(dates = dates,
                            dna = x$dna,
                            ctd = ctd,
                            w_dens = x$w,
                            ctd_directed = FALSE)

    data <- outbreaker_data(data = data)
    
    ## check the number of contacts are correct
    expect_equal(2*nrow(ctd), sum(data$contacts))

    ## make sure directionality is carried through
    expect_false(data$ctd_directed)

    

    ## identify non-matching labels
    wrong_dna <- x$dna
    rownames(wrong_dna) <- 1:length(x$onset)

    expect_error(data <- outbreaker_data(dates = dates,
                                         dna = wrong_dna,
                                         ctd = ctd,
                                         w_dens = x$w,
                                         ctd_directed = TRUE),
                 "DNA sequence labels don't match case ids")
    

})



