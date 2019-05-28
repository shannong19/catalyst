## SKG
## May 23, 2019
## LOL where does the time go
## Minimize joint MSE for SIR data for beta, gamma, and N
## Very specific for ebola data


#' MSE for ebola data
#' 
#' @param pars c(beta, gamma, N)
#' @param obs data frame with t, S, I, and R (these are the RAW numbers, we will add additional N)
#' @param t0 initial time start
#' @param T final time
#' @param prob_type 1 for KM 0 for RF
#' @param N default is NULL, meaning we will estimate N.  Otherwise, we use the N provided
#' return mean sum of square errors (averaged over time points)
mse_bgn <- function(pars, obs, t0 = 0, T = 100, prob_type = 1, N = NULL){
    ## Get beta, gamma, N
    beta <- pars[1]
    gamma <- pars[2]
    if(is.null(N)) N <- round(pars[3] *  10^5)
    ## Subset the observations to the proper time
    sub_inds <- which(obs$t >= t0 & obs$t < T)
    sub_obs <- obs[sub_inds,]
    N0 <- sum(sub_obs[1,])  ## This is from the raw data, pre-adjustment
    N1 <- max(0, N - N0)[1]
    sub_obs[,1] <- sub_obs[,1] + N1
    ## Get new init values
    x0 <- as.integer(sub_obs[1, c("S", "I", "R")])

    exp_X <- sirLoopGroups(x0, beta, gamma, T - t0, prob_type)
    obs_sir <- as.matrix(sub_obs[, c("S", "I", "R")])
    obs_sir_scale <- rescale_SIR(obs_sir, obs_sir)
    exp_X_scale <- rescale_SIR(exp_X, obs_sir)
    
    mse <- 1 / (T-t0) * sum((exp_X_scale - obs_sir_scale)^2)
    return(mse)
    


}



#' Rescale S, I, and R values between 0 and 1 (for observed data), predicted values can be higher or lower
#'
#' @param obs matrix of S, I, R
#' @param raw matrix of S, I, and R.  Will scale to these
#' @return rescaled observation matrix
rescale_SIR <- function(obs, raw){

    out <- do.call('cbind',
                   lapply(1:ncol(obs),
                          function(ii){
                              min_raw <- min(raw[,ii], na.rm = TRUE)
                              max_raw <- max(raw[,ii], na.rm = TRUE)
                              (obs[,ii] - min_raw) / (max_raw - min_raw)
                          }))
    return(out)

}

#' Calculate R0 from delta method
#'
#' @param mu vector (beta, gamma)
#' @param variance matrix of beta, gamma
#' @param T total time steps
#' @return variance of r0
r0_var <- function(mu, sigma, T){
    ##
    beta <- mu[1]
    gamma <- mu[2]
    sbb <- sigma[1,1]
    sgg <- sigma[2,2]
    sbg <- sigma[1,2]
    ##
    var <- 1 / T * (sbb/gamma^2 - sbg * (1 + beta) / gamma^3 + beta * sgg / gamma^4)
    return(var)

}
