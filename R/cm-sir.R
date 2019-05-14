## SKG
## CM estimates for the basic SIR
## May 14, 2019
## Instead of using the slow AM fun times
## Run a bunch of separate SIRs 


#' Simulate (many) CM-SIR model(s)
#'
#' @param L total number of simulations
#' @param total amount of time (0:(T-1))
#' @param X0_vec vector of size Gx3 (S1, I1, R1, S2, I2, R2, ..., SG, IG, RG)
#' @param par_vec vector of Gx2 (beta1, gamma1, ..., betaG, gammaG)
#' @param times_vec the time splits
#' @param prob_fxn either "KM" for Kermack-McKendrick or "RF" for Reed Frost
#' @param use_exp_X logical.  Whether to use random or non-random prob in the prob_fxn.  Default is TRUE
#' @param do_par logical.  Should we run in parallel?  Default is FALSE.  TBD
#' @param t0 in X data frame return values t0:(T-1)
#' @param interaction_type "heterog" for separate SIRs and "homog" for same
#' @return list with X_mat a  matrix where the columns are t, S, I, R, and ll
cm_sir <- function(L, T,
                   X0_vec,
                   par_vec,
                   times_vec,
                   prob_fxn = "KM",
                   use_exp_X = TRUE,
                   do_par = FALSE,
                   t0 = 0,
                   interaction_type = "heterog"){

    if(interaction_type == "heterog"){
        G <- length(times_vec)
        sims_X <- vector(mode = "list", length = L)
        N <- sum(X0_vec) / 3
        for(ll in 1:L){
            X_list <- vector(mode = "list", length = G)
            for(gg in 1:G){

                X0 <- X0_vec[((gg-1)  * 3 + 1):((gg-1) * 3 + 3)]
                N <- sum(X0)
                beta <- par_vec[(gg -1) * 2 + 1]
                gamma <- par_vec[(gg -1) * 2 + 2]
                if(gg == 1){
                    tstar <- t0
                } else {
                    tstar <- max(c(times_vec[gg-1], t0))
                }
                inf_prob <- get_CM_prob_inf(X0, beta, gamma,  prob_fxn,
                                            use_exp_X, T - tstar)

                X_mat <- update_X(X0, beta, gamma, inf_prob, tstar, T, ll)
                X_list[[gg]] <- X_mat
            }

            X <- combine_sims_X(X_list, N)
            sims_X[[ll]] <- X
        }
        out <- do.call('rbind', sims_X)
        out_list <- list(X = out)

    } else if(interaction_type == "homog"){

    } else{
        stop("Choose either 'heterog' or 'homog'")
    }
    return(out_list)

}
                   

#' Get the infection probabilities for the SIR
#' @param X0 vector (S0, I0, R0)
#' @param beta numeric
#' @param prob_fxn either "KM" or "RF"
#' @param use_exp_X logical.  If TRUE then we use I from SIRLoop.  Else we return NULL
#' @param T total time steps
#' @return vector of size T-1 which the probability of becoming infectious from time t-1 to t
get_CM_prob_inf <- function(X0, beta, gamma, prob_fxn,
                            use_exp_X, T){

    N <- sum(X0)
    inf_prob <- numeric(T-1)
    if(use_exp_X){
        pf <- ifelse(prob_fxn == "KM", 0, 1)
        X_out <- sirLoopGroups(X0, beta, gamma, T, pf)
        if(prob_fxn == "KM"){
            inf_prob <- beta * X_out[1:(T-1), 2] / N
        } else if(prob_fxn == "RF"){
            inf_prob <- ( 1- beta / N)^X_out[1:(T-1),2]
        }
    } else {
        inf_prob <- NULL
    }

    return(inf_prob)

}



#' update X for the CM SIR
#' @param X0 (S0, I0, R0)
#' @param beta numeric
#' @param gamma numeric
#' @param inf_prob either Null or vector of T - t0 - 1
#' @param t0 time we start at
#' @param T T-1 is the last time
#' @param ll simulation number
#' @return X matrix with columns t, S, I, R, ll nrow is T-1
update_X <- function(X0, beta, gamma, inf_prob, t0, T, ll){

    X <- matrix(0, nrow = T - t0, ncol = 5)
    N <- sum(X0)
    X[1,] <- c(t0, X0, ll)
    colnames(X) <- c("t", "S", "I", "R", "ll")
    for(tt in 1:(T-t0 - 1)){
        if(is.null(inf_prob)){
            ## Use actual X values to get inf_prob
            stop("Fill in")
        } else {
            cur_prob_inf <- inf_prob[tt]
        }
        X[tt + 1, 1] <- tt + t0 #t
        X[tt + 1, 2] <- X[tt, 2] - rbinom(1, X[tt, 2], prob = cur_prob_inf) 
        X[tt + 1, 4] <- X[tt, 4] + rbinom(1, X[tt, 3], prob = gamma)
        X[tt + 1, 3] <- N - X[tt + 1, 2] - X[tt + 1, 4]
    }
    X[, 5] <- ll
    return(X)

}
