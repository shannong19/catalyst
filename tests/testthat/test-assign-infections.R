## SKG
context("Testing assign_inf_recov_dates functions")

test_that("assign.byage", {
    
    full_df <- data.frame(AGE = c(2, 40, 12, 85, 85),
                          banana = c(T, T, T, F, T))
    inf_df <- data.frame(AGE = c(39, 2, 13, 85, 85),
                         I = c(1, 13, 2, 5, 20),
                         R = c(9, 15, 5, 14, 22))
    arg_list <- list(age_cutoff = 1, sequential = FALSE)

    out <- assign.byage(full_df = full_df,
                        inf_df = inf_df,
                        arg_list = arg_list)

    expect_equal(dim(out), c(5, 4))
    
    ## Test where there aren't enough matches
    full_df <- data.frame(AGE = c(2, 40, 12, 85, 19),
                          banana = c(T, T, T, F, T))
    inf_df <- data.frame(AGE = c(39, 2, 13, 85, 85),
                         I = c(1, 13, 2, 5, 20),
                         R = c(9, 15, 5, 14, 22))
    arg_list <- list(age_cutoff = 1, sequential = FALSE)
    out <- assign.byage(full_df = full_df,
                        inf_df = inf_df,
                        arg_list = arg_list)

    expect_equal(dim(out), c(5, 4))
})
