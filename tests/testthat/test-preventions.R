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


test_that("am_sir_one_sim with preventions",{
    ll <- 1
    T <- 5
    N <- 5
    A <- matrix(0, nrow = T, ncol = N)
    A[1, ] <- c(1, rep(0, N-1))
    beta <- .8
    gamma <- .2
    par1_vec <- rep(beta, N)
    gamma_vec <- rep(gamma, N)
    nbr_list <- NULL
    X_theor <- NULL

    out <- am_sir_one_sim(ll = ll, 
                          A = A,
                          par1_vec = par1_vec,
                          gamma_vec = gamma_vec,
                          nbr_list = nbr_list,
                          X_theor = X_theor)


    set.seed(2019)
    hh_mat <- matrix(sample(1:3, N, replace = TRUE), ncol = 1)
    preventions_list <- initializeNeighbors(hh_mat)
    preventions_type <- "isolation"
    env_df <- as.data.frame(hh_mat)
    colnames(env_df) <- c("hh")
    preventions_vars <- c("hh")

      out <- am_sir_one_sim(ll = ll, 
                          A = A,
                          par1_vec = par1_vec,
                          gamma_vec = gamma_vec,
                          nbr_list = nbr_list,
                          X_theor = X_theor,
                          do_preventions = TRUE,
                          preventions_type = preventions_type,
                          env_df = env_df,
                          preventions_vars = preventions_vars)

    
}) 

          
test_that("removeNeighbors", {
    nbr_list <- list('0'=-1, '1'=2, '2'=1, '3'=-1, '4'=-1)
    prev_list <- list('0'=-1, '1'=2, '2'=1, '3'=-1, '4'=-1)
    prev_inds <- 0:4
    inf_inds <- 0:4
    out_list <- removeNeighbors(nbr_list, prev_list, inf_inds, prev_inds)
    exp_list <- list('0'=-1, '1'=2, '2'=1, '3'=-1, '4'=-1)
    expect_equal(out_list, exp_list)
    ##
    nbr_list <- list('0'=-1, '1'=c(2,3), '2'=1, '3'=1, '4'=-1)
    prev_list <- list('0'=-1, '1'=3, '2'=-1, '3'=1, '4'=-1)
    prev_inds <- 0:4
    inf_inds <- 0:4
    out_list <- removeNeighbors(nbr_list, prev_list, inf_inds, prev_inds)
    exp_list <- list('0'=-1, '1'=3, '2'=-1, '3'=1, '4'=-1)
    expect_equal(out_list, prev_list)
    ##
    nbr_list <- list('0'=-1, '1'=c(2,3), '2'=1, '3'=1, '4'=-1)
    prev_list <- list('0'=-1, '1'=3, '2'=-1, '3'=1, '4'=-1)
    prev_inds <- 1
    inf_inds <- 0:4
    out_list <- removeNeighbors(nbr_list, prev_list, inf_inds, prev_inds)
    exp_list <- list('0'=-1, '1'=3, '2'=-1, '3'=1, '4'=-1)
    expect_equal(out_list, exp_list)
     ##
    nbr_list <- list('0'=-1, '1'=c(2,3), '2'=1, '3'=1, '4'=-1)
    prev_list <- list('0'=-1, '1'=c(2,3), '2'=1, '3'=1, '4'=-1)
    prev_inds <- 1
    inf_inds <- 0:4
    out_list <- removeNeighbors(nbr_list, prev_list, inf_inds, prev_inds)
    exp_list <- list('0'=-1, '1'=3, '2'=-1, '3'=1, '4'=-1)
    expect_equal(out_list, prev_list)


   
    })


test_that("update_preventions", {
    nbr_list <- list('0'=-1, '1'=2, '2'=1, '3'=-1, '4'=-1)
    prev_list <- list('0'=-1, '1'=2, '2'=1, '3'=-1, '4'=-1)
    prev_inds <- 0:4
    inf_inds <- 0:4
    days_infectious <- NULL
    preventions_delay <- 0
    out <- update_preventions(nbr_list,
                              prev_list,
                              prev_inds,
                              inf_inds,
                              days_infectious,
                              preventions_delay)
    expect_equal(out$new_nbrs, prev_list)
    expect_true(out$new_do_preventions)
    ##
    nbr_list <- list('0'=-1, '1'=2, '2'=1, '3'=-1, '4'=-1)
    prev_list <- list('0'=-1, '1'=2, '2'=1, '3'=-1, '4'=-1)
    prev_inds <- 0:4
    inf_inds <- 0:4
    days_infectious <- rep(0, 5)
    preventions_delay <- 0
    out <- update_preventions(nbr_list,
                              prev_list,
                              prev_inds,
                              inf_inds,
                              days_infectious,
                              preventions_delay)
    expect_equal(out$new_nbrs, prev_list)
    expect_true(out$new_do_preventions)
    ##
    nbr_list <- list('0'=-1, '1'=2, '2'=1, '3'=-1, '4'=-1)
    prev_list <- list('0'=-1, '1'=2, '2'=1, '3'=-1, '4'=-1)
    prev_inds <- 0:4
    inf_inds <- 2
    days_infectious <- rep(0, 5)
    preventions_delay <- 0
    out <- update_preventions(nbr_list,
                              prev_list,
                              prev_inds,
                              inf_inds,
                              days_infectious,
                              preventions_delay)
    expect_equal(out$new_nbrs, prev_list)
    expect_true(out$new_do_preventions)
        ##
    nbr_list <- list('0'=-1, '1'=2, '2'=1, '3'=-1, '4'=-1)
    prev_list <- list('0'=-1, '1'=2, '2'=1, '3'=-1, '4'=-1)
    prev_inds <- 0:4
    inf_inds <- 2
    days_infectious <- rep(0, 5)
    preventions_delay <- 1
    out <- update_preventions(nbr_list,
                              prev_list,
                              prev_inds,
                              inf_inds,
                              days_infectious,
                              preventions_delay)
    expect_equal(out$new_nbrs, nbr_list)
    expect_true(out$new_do_preventions)

    })
