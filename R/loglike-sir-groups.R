## SKG
## March 28, 2019
## Log likelihood functions for when we have many groups, blah



#' Calculate the loglike for data of SIR
#' 
#' @param par  vector of length 2xG where G is the number of groups.  The first G are beta_1, \dots, beta_g.  The second G are gamma_1, \dots, gamma_G
#' @param T number of steps 0:(T-1) inclusive
#' @param suff_stat either U a 3xN matrix with (A0, SMax, IMax) or A a TxN matrix where each entry is 0, 1, or 2 for S, I, and R, respectively
#' @param suff_stat_type character either "U" or "X" or "UX" default is U
#' @param prob_fxn character either "KM" for Kermack and McKendrick prob function or "RF" for Reed-Frost prob
#' @param neg_loglike logical.  Default is TRUE.  Should we return the negative loglike?
#' @param use_exp_X logical.  Default is TRUE.  Should we use expected values of X for the prob function or the data?
#' @param x0 c(S1, \dots, S_G, I1, \dots, I_G, R_1, \dots, R_G) , length Gx3
#' @param inf_nbrs TBD
#' @param X Tx(3G) matrix for use with loglike_sir.UX.  Default is NULL
#' @param Usub 3xN matrix, a partial U, for use with loglike_sir.UX.  Default is NULL
#' @param G_id vector of length N giving the group ID of each agent.  To be used with "U" or "UX" suff stats
#' @return scalar loglike
loglike_sirG <- function(par, T, suff_stat, suff_stat_type = "U",
                        prob_fxn = "KM", neg_loglike = TRUE,
                        use_exp_X = TRUE, x0 = NULL,
                        inf_nbrs = NULL,
                        X = NULL,
                        Usub = NULL){
   
    if(suff_stat_type == "U"){
        stop("Feature coming later")
        loglike <- loglike_sirG.U(par, T, suff_stat, prob_fxn,
                                 use_exp_X, x0,
                                 inf_nbrs,
                                 G_id)
    } else if(suff_stat_type == "X"){
        stop("Feature coming later")
        loglike <- loglike_sirG.X(par, T, suff_stat, prob_fxn,
                                 use_exp_X, x0,
                                 inf_nbrs)
    } else if(suff_stat_type == "UX"){
        loglike <- loglike_sirG.UX(par, T,
                                  Usub, X,
                                  prob_fxn = prob_fxn,
                                  use_exp_X = use_exp_X, x0 = X[1,],
                                  inf_nbrs = NULL,
                                  G_id)
        
    } else{
        stop("select a proper suff stat")
    }
    if(neg_loglike){
        loglike <- loglike * -1
    }

    return(sum(loglike))

}


#' @param par  vector of length 2xG where G is the number of groups.  The first G are beta_1, \dots, beta_g.  The second G are gamma_1, \dots, gamma_G
#' @param T number of steps 0:(T-1) inclusive
#' @param suff_stat either U a 3xN matrix with (A0, SMax, IMax) or A a TxN matrix where each entry is 0, 1, or 2 for S, I, and R, respectively
#' @param suff_stat_type character either "U" or "X" or "UX" default is U
#' @param prob_fxn character either "KM" for Kermack and McKendrick prob function or "RF" for Reed-Frost prob
#' @param neg_loglike logical.  Default is TRUE.  Should we return the negative loglike?
#' @param use_exp_X logical.  Default is TRUE.  Should we use expected values of X for the prob function or the data?
#' @param x0 c(S1, \dots, S_G, I1, \dots, I_G, R_1, \dots, R_G) , length Gx3
#' @param inf_nbrs TBD
#' @param X Tx(3G) matrix for use with loglike_sir.UX.  Default is NULL
#' @param Usub 3xN matrix, a partial U, for use with loglike_sir.UX.  Default is NULL
#' @param G_id vector of length N giving the group ID of each agent.  To be used with "U" or "UX" suff stats
#' @return scalar loglike
loglike_sirG.X <- function(par, T, suff_stat, prob_fxn,
                           use_exp_X, x0,
                           inf_nbrs){


}
