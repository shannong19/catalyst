## SKG
## March 6, 2019
##
## SIR Loglike functions
## 1. AM-SIR with U and KM
## 2. AM-SIR with U and RF
## 3. CM-SIR with X and KM
## 4. CM-SIR with X and RF

#' Calculate the loglike for data of SIR
#' 
#' @param par c("beta", "gamma")
#' @param T number of steps 0:(T-1) inclusive
#' @param suff_stat either U a 3xN matrix with (A0, SMax, IMax) or A a TxN matrix where each entry is 0, 1, or 2 for S, I, and R, respectively
#' @param suff_stat_type character either "U" or "A" default is U
#' @param prob_fxn character either "KM" for Kermack and McKendrick prob function or "RF" for Reed-Frost prob
#' @param neg_loglike logical.  Default is TRUE.  Should we return the negative loglike?
#' @param use_exp_X logical.  Default is TRUE.  Should we use expected values of X for the prob function or the data?
#' @param x0 c(S0, I0, R0)
#' @param inf_nbrs TBD
#' @return scalar loglike
loglike_sir <- function(par, T, suff_stat, suff_stat_type = "U",
                        prob_fxn = "KM", neg_loglike = TRUE,
                        use_exp_X = TRUE, x0 = NULL,
                        inf_nbrs = NULL){
   
    if(suff_stat_type == "U"){
        loglike <- loglike_sir.U(par, T, suff_stat, prob_fxn,
                                 use_exp_X, x0,
                                 inf_nbrs)
    } else if(suff_stat_type == "A"){
        loglike <- loglike_sir.A(par, T, suff_stat, prob_fxn,
                                 use_exp_X, x0,
                                 inf_nbrs)
    } else{
        stop("select a proper suff stat")
    }
    if(neg_loglike){
        loglike <- loglike * -1
    }

    return(sum(loglike))

}


#' Get loglike for SIR with U statistics
#' 
#' @param par c("beta", "gamma")
#' @param T number of steps 0:(T-1) inclusive
#' @param U  U a 3xN matrix with (A0, SMax, IMax) 
#' @param suff_stat_type character either "U" or "A" default is U
#' @param fxn function either "KM" for Kermack and McKendrick prob function or "RF" for Reed-Frost prob
#' @param neg_loglike logical.  Default is TRUE.  Should we return the negative loglike?
#' @param use_exp_X logical.  Default is FALSE.  Should we use expected values of X for the prob function or the data?
#' @param x0 c(S0, I0, R0) or NULL
#' @param inf_nbrs TBD
#' @return scalar loglike
loglike_sir.U <- function(par, T, U, prob_fxn = "KM",
                          use_exp_X = TRUE, x0 = NULL,
                          inf_nbrs = NULL){
    if(prob_fxn == "KM"){
        fxn <- KM_prob
    } else if(prob_fxn == "RF"){
        fxn <- RF_prob
    } else{
        stop("select a proper prob_fxn")
    }
    ##
    if(use_exp_X){
        pf <- ifelse(prob_fxn == "KM", 0, 1)
        X <- sirLoop(x0,
                     par[1], par[2], T, pf)
    } else{
        X <- UtoX_SIR(U, T)
    }
    N <- sum(X[1,])
    pt <- fxn(par, X, inf_nbrs, N)
    ## Calc loglike for each agent
    loglike <- sapply(1:N, function(ii){
        loglike_agent_sir(pt, par[2], T, U[,ii])
    })
    return(loglike)

}

#' Get loglike for an individual agent
#'
#' @param pt vector of size T-1.  entry i is probablity of becoming infectious from time i-1 to i for i=1, \dots, T-1.  Remember time starts at 0
#' @param gamma probability of recovery
#' @param T max time steps
#' @param Un suff stat for agent (A0, SMax, IMax).  Remember SMax, IMax start at 0.
#' @return loglike for agent
loglike_agent_sir <- function(pt, gamma, T, Un){
    A0 <- Un[1]
    SMax <- Un[2] + 1
    IMax <- Un[3] + 1
    pSMax <- ifelse(SMax == T, 1, pt[SMax]) ## all to avoid an out of bounds issue
    if(A0 == 0){
        loglike <- (SMax > 1) * sum(log( 1 - pt[1:(SMax-1)])) +
            (SMax < T) * log(pSMax) +
            (SMax < T) * (IMax - SMax - 1) * log( 1- gamma) +
            (IMax == T & SMax < (T-1)) * log(1 - gamma) +
            (IMax < T) * log(gamma)
    } else if(A0 == 1){
        loglike <- (IMax - 1) * log(1 - gamma) +
            (IMax < T) * log(gamma)
    } else if(A0 == 2){
        loglike <- 0
    } else{
        stop("This is not an applicable state")
    }
    return(loglike)

}
