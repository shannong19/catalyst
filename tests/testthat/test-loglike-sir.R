## SKG
context("loglike-sir")

test_that("loglike_agent sir", {

    T <- 4
    pt <- seq(.25, .75, length.out = 4)
    gamma <- .2
    Un <- c(0, 0, 1) # is iniitally S, becomes I from time 0 to 1.  Becomes R from 1 to 2
    out <- loglike_agent_sir(pt, gamma, T, Un)
    exp_out <- log(pt[1]) + log(gamma)
    expect_equal(out, exp_out)
    ##
    Un <- c(0, 1, 3) # is iniitally S, becomes I from time 1 to 2.  Becomes R never
    out <- loglike_agent_sir(pt, gamma, T, Un)
    exp_out <- sum(log(1-pt[1]))  + log(pt[2]) +  2 * log(1 -gamma)
    expect_equal(out, exp_out)
    ##
    Un <- c(0, 2, 3) # is iniitally S, becomes I from time 2 to 3.  Becomes R never
    out <- loglike_agent_sir(pt, gamma, T, Un)
    exp_out <- sum(log(1-pt[1:2]))  + log(pt[3]) 
    expect_equal(out, exp_out)
    ##
    Un <- c(0, 3, 3) # is iniitally S, becomes I never.  Becomes R never
    out <- loglike_agent_sir(pt, gamma, T, Un)
    exp_out <- sum(log(1-pt[1:3])) 
    expect_equal(out, exp_out)
    ##
    Un <- c(0, 0, 2) # is iniitally S, becomes I from 0 to 1.  Becomes R from 2 to 3
    out <- loglike_agent_sir(pt, gamma, T, Un)
    exp_out <- log(pt[1]) + (log(1 - gamma))  + log(gamma)
    expect_equal(out, exp_out)
    ##
    Un <- c(1, 0, 2) # is iniitally I, becomes I from 0 to 1.  Becomes R from 2 to 3
    out <- loglike_agent_sir(pt, gamma, T, Un)
    exp_out <-  2 * (log(1 - gamma))  + log(gamma)
    expect_equal(out, exp_out)
    ##
    Un <- c(1, 2, 2) # is iniitally I, becomes I from 0 to 1.  Becomes R from 2 to 3
    out <- loglike_agent_sir(pt, gamma, T, Un)
    exp_out <-  2 * (log(1 - gamma))  + log(gamma)
    expect_equal(out, exp_out)
    ##
    Un <- c(1, 2, 3) # is iniitally I, becomes I from 0 to 1.  Becomes R from 2 to 3
    out <- loglike_agent_sir(pt, gamma, T, Un)
    exp_out <-  3 * (log(1 - gamma))  
    expect_equal(out, exp_out)
    ##
    Un <- c(2, 2, 2) # is iniitally I, becomes I from 0 to 1.  Becomes R from 2 to 3
    out <- loglike_agent_sir(pt, gamma, T, Un)
    exp_out <-0
    expect_equal(out, exp_out)



})


test_that("loglike_sir.U", {
    U <- matrix(c(0, 0, 1, 0,
                  0, 1, 2, 2,
                  1, 2, 2, 2), byrow = TRUE, nrow = 3)
    T <- 3
    par <- c(.2, .7)
    X <- UtoX_SIR(U, T)
    x <- X[1,]
    out <- loglike_sir.U(par, T, U,
                         prob_fxn = "KM",
                         use_exp_X = TRUE,
                         x0 = x)
    expect_true(all(!is.na(out)))

})


test_that("loglike_sir", {
    par <- c(.2, .7)
    T <- 3
    U <- matrix(c(0, 0, 1, 0,
                  0, 1, 2, 2,
                  1, 2, 2, 2), byrow = TRUE, nrow = 3)
    X <- UtoX_SIR(U, T)
    x <- X[1,]
    out <- loglike_sir(par = par, T = T,
                       suff_stat = U, suff_stat_type = "U",
                       prob_fxn = "RF",
                       neg_loglike = TRUE,
                       use_exp_X = TRUE,
                       x0 = x,
                       inf_nbrs = NULL)
    expect_true(out > 0)
    

})

test_that("loglike_sir", {
    x <- c(950, 50, 0)
    par <- c(.1, .03)
    T <- 100
    X <- sirLoop(x, par[1], par[2], 100, 0)
    N <- sum(X[1,])

    plot(1, 1, xlim =c(1,T), ylim = c(0,N), type = "n")
    lines(1:T, X[,1], type = "l", lwd = 2, col = "blue")
    lines(1:T, X[,2], lwd = 2, col = "red")
    lines(1:T, X[,3], lwd = 2, col = "yellow3")

    
    out <- loglike_sir(par = par, T = T,
                       suff_stat = U, suff_stat_type = "U",
                       prob_fxn = "KM",
                       neg_loglike = TRUE,
                       use_exp_X = TRUE,
                       x0 = x0,
                       inf_nbrs = NULL)
   
    


})


test_that("sir", {
      
    L <- 10
    T <- 100
    N <-  1000
    A0 <- c(rep(0, .95 * N), rep(1, .05 * N))
    x0 <- c(950, 50, 0)
    prob_fxn <- "KM"
    par1_vec <- rep(.5, N)
    gamma_vec <- rep(0.25, N)
    ## par1_vec <- runif(N)
    ## gamma_vec <- runif(N)
    nbr_list <- NULL
    use_exp_X <- TRUE
    keep_A <- FALSE
    keep_U <- TRUE
    write_sim <- FALSE
    writing_list <- NULL
    do_par <- FALSE
    ##
    t <- proc.time()[3]
    sims <- am_sir(L = L, T = T,
                   A0 = A0, prob_fxn = prob_fxn,
                   par1_vec = par1_vec,
                   gamma_vec = gamma_vec,
                   nbr_list = nbr_list,
                   use_exp_X = use_exp_X,
                   keep_A = keep_A,
                   keep_U = keep_U,
                   write_sim = write_sim,
                   writing_list = NULL,
                   do_par = do_par)
    proc.time()[3] - t

    U <- sims$U[3,,]
    dim(U)        
    out <- optim(par = c(.2, .1), loglike_sir, T = T,
                 suff_stat = U, suff_stat_type = "U",
                 prob_fxn = "KM",
                 neg_loglike = TRUE,
                 use_exp_X = TRUE,
                 x0 = x0,
                 inf_nbrs = NULL)
    out$par

})


test_that("loglike_sir.X and friends", {


    par <- c(.2, .7)
    T <- 3
    pt <- rep(par[1], T-1)
    U <- matrix(c(0, 0, 1, 0,
                  0, 1, 2, 2,
                  1, 2, 2, 2), byrow = TRUE, nrow = 3)
    X <- UtoX_SIR(U, T)
    x <- X[1,]
    N <- sum(x)


    out <- loglike_SIR_CM(T, pt, par[2], X, N)
    exp_loglike <- (1 * log(pt[1]) + 2 * log(1-pt[1])) +
        (0 + 1 * log(1-par[2])) +
        (1 * log(pt[2]) + 1 * log(1-pt[2])) +
        (1 * log(par[2]) + 1 * log(1-par[2]))
    expect_equal(out, exp_loglike)


    
    out <- loglike_sir(par = par, T = T,
                       suff_stat = X, suff_stat_type = "X",
                       prob_fxn = "KM",
                       neg_loglike = TRUE,
                       use_exp_X = TRUE,
                       x0 = x,
                       inf_nbrs = NULL)

    
})


test_that("loglike_sir.UX and friends", {
    ## This is for when X is total and U is partial
    ## e.g. different agents in U are associated with different parameters, but still encounter the total number of infections in X

    par <- c(.2, .7)
    T <- 3
    U <- matrix(c(0, 0, 1, 0,
                  0, 1, 2, 2,
                  1, 2, 2, 2), byrow = TRUE, nrow = 3)
    X <- UtoX_SIR(U, T)
    x <- X[1,]
    Usub <- U[, 1:3]
    out2 <- loglike_sir.UX(par = par, T = T,
                          X = X, U= Usub,
                          prob_fxn = "KM",
                          use_exp_X = FALSE,
                          x0 = x,
                          inf_nbrs = NULL)
    expect_true(all(out2 < 0))

    out2 <- loglike_sir(par = par, T = T,
                       suff_stat = NULL, suff_stat_type = "UX",
                       X = X, Usub= Usub,
                       prob_fxn = "KM",
                       neg_loglike = TRUE,
                       use_exp_X = FALSE,
                       x0 = x,
                       inf_nbrs = NULL)
    expect_true(all(out2 > 0))

    ########################################
### Trying with multiple groups
###################################3

    G <- 2
    par <- c(.2, 1, 0, .7)
    G_id <- c(1, 1, 2, 1)
    T <- 3
    U <- matrix(c(0, 0, 1, 0,
                  0, 1, 2, 2,
                  1, 2, 2, 2), byrow = TRUE, nrow = 3)
    X <- UtoX_SIR(U, T)
    x <- X[1,]
    Usub <- U
    out <- loglike_sir.UX(par = par, T = T,
                          X = X, U= U,
                          prob_fxn = "KM",
                          use_exp_X = TRUE,
                          x0 = x,
                          inf_nbrs = NULL,
                          G_id = G_id,
                          G = G)

    

})


test_that("sir_init_groups", {
    A0 <- c(0, 0, 0, 1)
    G_id <- c(1, 2, 2, 1)
    G <- 2
    out <- sir_init_groups(A0, G_id, G)
    out
    expect_equal(out, c(1, 2, 1, 0, 0, 0))

    })



test_that("aggregate_X", {
    X <- matrix(c(9, 7, 3, 1,
                  1, 2, 5, 3,
                  0, 1, 2, 6), ncol = 3)
    out <- aggregate_X(X)
    expect_equal(sum(out - X), 0)
    ## Many groups
    X <- matrix(c(9, 7, 3, 1,
                  1, 2, 5, 3,
                  0, 1, 2, 6,
                  9, 7, 3, 1,
                  1, 2, 5, 3,
                  0, 1, 2, 6), ncol = 6)
    out <- aggregate_X(X)
    expect_equal(dim(out), c(4, 3))
    
    })
