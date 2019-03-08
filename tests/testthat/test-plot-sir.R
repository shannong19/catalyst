context("time to plot")

test_that("plot_sims", {
 
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



    g <- plot_sims.reg(sims, K = 3,
                   beta = .1,
                   gamma = .03,
                   N = N,
                   S0 = 950,
                   I0 = 50,
                   L = L,
                   T = 100,
                   pretty = TRUE)

    ##devtools::load_all()
    g <- plot_mean_var.reg(sims, K = 3,
                   beta = .1,
                   gamma = .03,
                   N = N,
                   S0 = 950,
                   I0 = 50,
                   L = L,
                   T = 100,
                   pretty = TRUE)


    ## Loglinear
    g <- plot_sims.loglin(sims, K = 3,
                   beta = .1,
                   gamma = .03,
                   N = N,
                   S0 = 950,
                   I0 = 50,
                   L = L,
                   T = 100,
                   pretty = TRUE)
    

})


test_that("plot_ests.state", {
    
    library(ggplot2)
    library(reshape2)
    
    t <- rep(1:100, 1)
    S <- 3 * t
    I <- log(t)
    R <- 2 + t
    df1 <- data.frame(t = t, S = S, I = I, R = R)
    dfm1 <- melt(df1, id.vars = "t")
    colnames(dfm1) <- c("t", "state", "obs")
    dfm1$mean = NA
    dfm1$var = NA
    dfm1$data_type = 0
    var_order <- c("t", "obs", "mean", "var", "data_type", "state")
    dfm1 <- dfm1[, var_order]
    ## est
    df2 <- data.frame(t = rep(t, 3),
                      mean = c(2.5 * S,
                               1.2 * log(I),
                               R,
                               2.2 * S,
                               1.4 * log(I),
                               R + 1),
                      var = c(rnorm(2 * 300, 300, 1)),
                              data_type = rep(c("1", "2"),
                                              each = 300 ),
                      state = rep(c("S", "I", "R",
                                    "S", "I", "R"), each = 100),
                      obs = NA)
    df2 <- df2[, var_order]

    df <- rbind(dfm1, df2)
    g <- plot_ests.state(df)
})
