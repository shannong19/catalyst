context("AM-SIR")

test_that("AM SIR functions",{
    
    ll <- 1
    prob_fxn <- "KM"
    nbr_list <- NULL
    X_theor <- NULL
    keep_A <- TRUE
    keep_U <- TRUE
    write_sim <- FALSE
    writing_list <- NULL
    N <- 4
    T <- 3
    A <- matrix(c(0, 0, 0, 1,
                  1, 0, 0, 1,
                  2, 0, 1, 1), byrow = TRUE, nrow = T, ncol = N)
    par1_vec <- rep(1,N)
    gamma_vec <- rep(0,N)

    ## t <- proc.time()[3]
    out <- am_sir_one_sim(ll = ll,
                      A = A,
                      prob_fxn = prob_fxn,
                      par1_vec = par1_vec,
                      gamma_vec = gamma_vec,
                      nbr_list = nbr_list,
                      X_theor= X_theor,
                      keep_A = keep_A,
                      keep_U = keep_U,
                      write_sim = write_sim,
                      writing_list = writing_list)
    ## proc.time()[3] - t
    expect_equal(length(out), 3)




    
})


test_that("full simulations", {
    
    L <- 100
    T <- 100
    N <-  1000
    A0 <- c(rep(0, .95 * N), rep(1, .05 * N))
    prob_fxn <- "KM"
    par1_vec <- rep(.1,N)
    gamma_vec <- rep(0.03,N)
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

    library(ggplot2)
    library(reshape2)
    x <- as.data.frame(sims$X)
    xm <- melt(x, id.vars = c("t", "ll"))
    head(xm)
    ggplot(xm, aes(x = t, y = value, group = paste(ll, variable),
                   col = variable)) +
        geom_line() + scale_color_manual(values = c("blue", "red", "yellow3"))
    
    

    
})
