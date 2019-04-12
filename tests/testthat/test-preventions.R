cotnext("Testing prevention functions")

test_that("precompute_preventions_nbrs", {

    N <- 5
    E <- 1
    env_df <- data.frame(HH_id = c(0, 1, 1, 3, 4))
    vars <- c("HH_id")

    out <- precompute_preventions_nbrs(env_df, vars)
    exp_out <- list('0'=-1, '1'=2, '2'=1, '3'=-1, '4'=-1)
    expect_equal(out, exp_out)

})


test_that("updateNeighbors", {
    nbr_list <- list('0'=-1, '1'=2, '2'=1, '3'=-1, '4'=-1)
    prev_list <- list('0'=-1, '1'=2, '2'=1, '3'=-1, '4'=-1)
    inds <- 0:4
    out_list <- updateNeighbors(nbr_list, prev_list, inds)
    expect_equal(out_list, nbr_list)
    ##
    prev_list <- list('0'=-1, '1'=c(2:4), '2'=c(1:4), '3'=-1, '4'=-1)
    out_list <- updateNeighbors(nbr_list, prev_list, inds)
    expect_equal(out_list, prev_list)
    ##
    inds <- c(0, 1, 3, 4)
    out_list <- updateNeighbors(nbr_list, prev_list, inds)
    expect_equal(out_list, list('0'=-1, '1'=c(2:4), '2'=1, '3'=-1, '4'=-1))

    mat <- matrix(rep(1, 5), ncol = 1)
    nbr_list <- initializeNeighbors(mat)
})
