context("Testing CM-SIR functions")

test_that("update X", {
    X0 <- c(9, 1, 0)
    beta <- .1
    gamma <- 0
    inf_prob <- c(0, 0, 1)
    T <- 4
    t0 <- 0
    ll <- 1
    out <- update_X(X0, beta, gamma, inf_prob, t0, T, ll)
    exp_out <- matrix(c(9, 1, 0,
                        9, 1, 0,
                        9, 1, 0,
                        0, 10, 0), byrow = TRUE, ncol = 3)
    colnames(exp_out) <- c("S", "I", "R")
    expect_equal(out[, 2:4], exp_out)
    


})


test_that("get_CM_prob_inf",{
     X0 <- c(950, 50, 0)
     beta <- .1
     gamma <- .03
     prob_fxn <- "KM"
     use_exp_X <- TRUE
     T <- 100
     prob_inf <- get_CM_prob_inf(X0, beta, gamma,
                                 prob_fxn, use_exp_X, T)
     expect_equal(length(prob_inf), T - 1)


     ###########################
     X0 <- c(9, 1, 0)
     beta <- 1
     gamma <- 0
     prob_fxn <- "KM"
     use_exp_X <- TRUE
     T <- 5
     prob_inf <- get_CM_prob_inf(X0, beta, gamma,
                                 prob_fxn, use_exp_X, T)
     exp_out <- .1
     expect_equal(prob_inf[1], exp_out)

     ##
     X0 <- c(9, 1, 0)
     beta <- 0
     gamma <- 1
     prob_fxn <- "KM"
     use_exp_X <- TRUE
     T <- 5
     prob_inf <- get_CM_prob_inf(X0, beta, gamma,
                                 prob_fxn, use_exp_X, T)
     exp_out <- rep(0, 5-1)
     expect_equal(prob_inf, exp_out)

     ##
     X0 <- c(9, 1, 0)
     beta <- 0
     gamma <- 1
     prob_fxn <- "KM"
     use_exp_X <- FALSE
     T <- 5
     prob_inf <- get_CM_prob_inf(X0, beta, gamma,
                                 prob_fxn, use_exp_X, T)
     expect_true(is.null(prob_inf))

     
})


test_that("cm_sir", {
    
    L <- 10
    X0_vec <- c(950, 50, 0)
    par_vec <- c(.1, .03)
    T <- 100
    times_vec <- T
    out <- cm_sir(L = L, T = T,
                  X0_vec = X0_vec,
                  par_vec = par_vec,
                  times_vec = times_vec,
                  prob_fxn = "KM",
                  use_exp_X = TRUE,
                  do_par = FALSE,
                  t0 = 0,
                  interaction_type = "heterog")
    expect_equal(dim(out[[1]]), c(T * L, 5))

    ## MULTIPLE groups
    L <- 100
    X0_vec <- c(950, 49, 0, 0, 1, 0)
    par_vec <- c(.1, .03, .25, .2)
    T <- 100
    times_vec <- c(10, T)
    out <- cm_sir(L = L, T = T,
                  X0_vec = X0_vec,
                  par_vec = par_vec,
                  times_vec = times_vec,
                  prob_fxn = "KM",
                  use_exp_X = TRUE,
                  do_par = FALSE,
                  t0 = 0,
                  interaction_type = "heterog")
    out2 <- am_plot_mean_var(X = out$X, obs = NULL,
                            plot_var = TRUE,
                            cols = c("black", "blue"),
                            title = "woo")



})
head(out2$plot_df)
