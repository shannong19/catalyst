context("SSE and MSE")

test_that("sse_sir.df", {
    T <- 4
    N <- 4
    obs <- data.frame(t = 1:T,
                      S = c(3, 2, 2, 1),
                      I = c(1, 1, 1, 2),
                      R = c(0, 1, 1, 1))
    exp <- data.frame(t = 1:T,
                      S = c(3, 2, 1, 1),
                      I = c(1, 2, 1, 0),
                      R = c(0, 0, 2, 3))

    out <- sse_sir.df(obs, exp)
    expect_equal(sum((obs - exp)^2), as.integer(out[1]))

    ##
    obs <- data.frame(t = c(1, 3, 4),
                      S = c(3, 2, 1),
                      I = c(1, 1, 2),
                      R = c(0, 1, 1))
    exp <- data.frame(t = 1:T,
                      S = c(3, 2, 1, 1),
                      I = c(1, 2, 1, 0),
                      R = c(0, 0, 2, 3))

    out <- sse_sir.df(obs, exp)
    tot <- 0 + 2 + 8 
    expect_equal(12 -2 , as.integer(out[1]))

    ##
        ##
    exp <- data.frame(t = c(1, 3, 4),
                      S = c(3, 2, 1),
                      I = c(1, 1, 2),
                      R = c(0, 1, 1))
    obs <- data.frame(t = 1:T,
                      S = c(3, 2, 1, 1),
                      I = c(1, 2, 1, 0),
                      R = c(0, 0, 2, 3))
    out <- sse_sir.df(obs, exp)
    tot <- 0 + 2 + 8 
    expect_equal(tot , as.integer(out[1]))

    ##
    ##
    exp <- data.frame(t = c(1, 3),
                      S = c(3, 2),
                      I = c(1, 1),
                      R = c(0, 1))
    obs <- data.frame(t = 2:T,
                      S = c(2, 1, 1),
                      I = c(2, 1, 0),
                      R = c(0, 2, 3))
    out <- sse_sir.df(obs, exp)
    tot <- 2 
    expect_equal(tot , as.integer(out[1]))
    expect_equal(1, as.integer(out[2]))


})


test_that("sse_sir", {
    exp <- data.frame(t = c(1, 2),
                      S = c(3, 2),
                      I = c(1, 1),
                      R = c(0, 1))
    obs <- data.frame(t = c(1, 2, 1, 2),
                      S = c(2, 1, 1, 0),
                      I = c(2, 1, 3, 1),
                      R = c(0, 2, 0, 3),
                      ll = c(1, 1, 2, 2))
    out <- sse_sir(obs, exp)
    tot <- (2 + 2) + (8 + 8)
    expect_equal(tot , as.integer(out[1]))
    expect_equal(4, as.integer(out[2]))
    ##
    exp <- data.frame(t = c(1, 2),
                      S = c(3, 2),
                      I = c(1, 1),
                      R = c(0, 1))
    obs <- data.frame(t = c(1, 2, 1, 3),
                      S = c(2, 1, 1, 0),
                      I = c(2, 1, 3, 1),
                      R = c(0, 2, 0, 3),
                      ll = c(1, 1, 2, 2))
    out <- sse_sir(obs, exp)
    tot <- (2 + 2) + (8 + 0)
    expect_equal(tot , as.integer(out[1]))
    expect_equal(3, as.integer(out[2]))
})


test_that("min_mse_obj", {
    x0 <- c(950, 50, 0)
    T <- 100
    prob_type <- 0 # KM
    par <- c(.1, .3)
    X <- sirLoop(x0, par[1], par[2], T, prob_type)
    df <- as.data.frame(X)
    colnames(df) <- c("S", "I", "R")
    df$t <- 0:(T-1)
    df$ll <- 1

    out <- min_mse_obj(par = par, data = df,
                       x0 = x0,
                       T = T,
                       prob_type = prob_type)
    expect_equal(out, 0)


    ##
    x0 <- c(950, 50, 0)
    T <- 50
    prob_type <- 0 # KM
    par <- c(.1, .3)
    X <- sirLoop(x0, par[1], par[2], T, prob_type)
    df <- as.data.frame(X)
    colnames(df) <- c("S", "I", "R")
    df$t <- 0:(T-1)
    df$ll <- 1

    out <- min_mse_obj(par = par, data = df,
                       x0 = x0,
                       T = 100,
                       prob_type = prob_type)
    expect_equal(out, 0)

})


test_that("full scale simulation", {
    
    L <- 10
    T <- 100
    N <-  1000
    A0 <- c(rep(0, .95 * N), rep(1, .05 * N))
    prob_fxn <- "KM"
    bet <- .1
    gamma <- .03
    par1_vec <- rep(beta, N)
    gamma_vec <- rep(gamma ,N)
    ## par1_vec <- runif(N)
    ## gamma_vec <- runif(N)
    nbr_list <- NULL
    use_exp_X <- TRUE
    keep_A <- FALSE
    keep_U <- FALSE
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
    

    x0 <- c(950, 50, 0)
    T <- 100
    prob_type <- 0
    ## Not fast
    out <- optim(par = c(.3, .1), fn = min_mse_obj, data = as.data.frame(sims$X),
                 x0 = x0, T = T, prob_type = 1)
    expect_equal(1, 1)
})
