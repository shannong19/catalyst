## SKG
## March 7, 2019
## MSE for the SIR.  duh


#' Calculate the joint sum of squared errors for S, I, and R
#'
#' @param obs data frame with columns t, S, I, R, and ll: the time, #S, #I, #R and simulation number, respectively
#' @param exp  data frame with columns t, S, I, R: the time, #S, #I, #R and simulation number, respectively
#' @return vector of sum of squared errors and effective number of observations (total number of observations where there is also a recorded point in the expected df)
#' @details expects that if 
sse_sir <- function(obs, exp){
    sse <- plyr::daply(obs, .var = c("ll"), .fun = sse_sir.df, exp = exp)
    return(colSums(matrix(sse, ncol = 2)))

}

#' Calculate the joint sum of squared errors for S, I, and R
#'
#' @param obs data frame with columns t, S, I, R, and ll: the time, #S, #I, #R and simulation number, respectively.  There is only one unique value for ll
#' @param exp  data frame with columns t, S, I, R: the time, #S, #I, #R and simulation number, respectively
#' @return sum of squared errors and number of obsrvations used to compute it
#' ##TODO: put in c++
sse_sir.df <- function(obs, exp){
    colnames(exp)[1:4] <- c("t", "S_exp", "I_exp", "R_exp")
    df_join <- plyr::join(obs, exp, by = "t")
    sse_S <- sum((df_join$S  - df_join$S_exp)^2, na.rm = TRUE)
    sse_I <- sum((df_join$I  - df_join$I_exp)^2, na.rm = TRUE)
    sse_R <- sum((df_join$R  - df_join$R_exp)^2, na.rm = TRUE)
    sse <- sse_S + sse_I + sse_R
    n_obs <- nrow(na.omit(df_join))
    return(c(sse = sse, n_obs = n_obs))

}


#' Objective function set up for finding best params
#' @param par (beta, gamma) or (rho, gamma)
#' @param data data frame with columns t, S, I, R, and ll: the time, #S, #I, #R and simulation number, respectively.  There is only one unique value for ll.  Assumes it will match up for time t=0 for the expected value
#' @param x0 vector of initial (S0, I0, R0)
#' @param T total time steps 0:(T-1)
#' @param prob_type either 0 or 1 for KM prob, or Reed Frost, respectively
#' @return joint MSE 
min_mse_obj <- function(par = c(.5, .25), data, x0, T,
                        prob_type = 0){
    X <- sirLoop(x0, par[1], par[2], T, prob_type)
    exp <- data.frame(t= 0:(T-1), X)
    colnames(exp) <- c("t", "S", "I", "R")
    out <- sse_sir(obs = data, exp = exp)
    return(as.numeric(out[1] / out[2]))


}
