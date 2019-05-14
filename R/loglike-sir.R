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
#' @param suff_stat_type character either "U" or "X" or "UX" default is U
#' @param prob_fxn character either "KM" for Kermack and McKendrick prob function or "RF" for Reed-Frost prob
#' @param neg_loglike logical.  Default is TRUE.  Should we return the negative loglike?
#' @param use_exp_X logical.  Default is TRUE.  Should we use expected values of X for the prob function or the data?
#' @param x0 c(S0, I0, R0)
#' @param inf_nbrs TBD
#' @param X Tx3 matrix for use with loglike_sir.UX.  Default is NULL
#' @param Usub 3xNstar matrix, a partial U, for use with loglike_sir.UX.  Default is NULL
#' @param G_id group ID.  Default is NULL
#' @param G total number of groups.  Default is 1
#' @return scalar loglike
loglike_sir <- function(par, T, suff_stat, suff_stat_type = "U",
                        prob_fxn = "KM", neg_loglike = TRUE,
                        use_exp_X = TRUE, x0 = NULL,
                        inf_nbrs = NULL,
                        X = NULL,
                        Usub = NULL,
                        G_id = NULL,
                        G = 1){
   
    if(suff_stat_type == "U"){
        loglike <- loglike_sir.U(par, T, suff_stat, prob_fxn,
                                 use_exp_X, x0,
                                 inf_nbrs)
    } else if(suff_stat_type == "X"){
        loglike <- loglike_sir.X(par, T, suff_stat, prob_fxn,
                                 use_exp_X, x0,
                                 inf_nbrs)
    } else if(suff_stat_type == "UX"){
        loglike <- loglike_sir.UX(par, T,
                                  Usub, X,
                                  prob_fxn = prob_fxn,
                                  use_exp_X = use_exp_X, x0 = X[1,],
                                  inf_nbrs = NULL,
                                  G_id = G_id,
                                  G = G)
        
    } else{
        stop("select a proper suff stat")
    }
    if(neg_loglike){
        loglike <- loglike * -1
    }

    return(sum(loglike))

}


#' Get loglike for SIR with partial U statistics and full X
#' 
#' @param par c("beta", "gamma")
#' @param T number of steps 0:(T-1) inclusive
#' @param U  U a 3xNstar matrix with (A0, SMax, IMax)
#' @param X Tx3 matrix
#' @param suff_stat_type character either "U" or "X" default is U.
#' @param fxn function either "KM" for Kermack and McKendrick prob function or "RF" for Reed-Frost prob
#' @param neg_loglike logical.  Default is TRUE.  Should we return the negative loglike?
#' @param use_exp_X logical.  Default is FALSE.  Should we use expected values of X for the prob function or the data?
#' @param x0 c(S0, I0, R0) or NULL
#' @param inf_nbrs TBD
#' @param G_id group ID.  Default is NULL
#' @param G total number of groups.  Default is 1
#' @return scalar loglike
loglike_sir.UX <- function(par, T,
                           U, X,
                           prob_fxn = "KM",
                           use_exp_X = FALSE, x0 = NULL,
                           inf_nbrs = NULL,
                           G_id = NULL,
                           G = 1){
    N <- sum(X[1,])
    betas <- par[1:G]
    gammas <- par[(G+1):(2 * G)]
    beta_n <- group_to_agent_par(betas, G_id, ncol(U))
    gamma_n <- group_to_agent_par(gammas, G_id, ncol(U))
    if(prob_fxn == "KM"){
        fxn <- KM_probG
        prob_type <- 0
    } else if(prob_fxn == "RF"){
        fxn <- RF_probG
        prob_type <- 1
    } else{
        stop("select a proper prob_fxn")
    }
    ##
    if(use_exp_X){
        x0 <- sir_init_groups(U[1,], G_id, G, X)  ## WARNING:  Assumes size of U is N, otherwise will give weird estimates
        X_group <- sirLoopGroups(x0,
                                 betas, gammas,
                                 T, prob_type)
        X <- aggregate_X(X_group)
    } else{
        X <- X
    }
    ## Constraining the pars
    beta_n <- ifelse(beta_n> 1, 1, beta_n)
    gamma_n<- ifelse(gamma_n<= 0, 1e-8, gamma_n)

    pt <- fxn(beta_n, X, inf_nbrs, N)
    ## Calc loglike for each agent
    loglike <- sapply(1:ncol(U), function(ii){
        loglike_agent_sir(pt[,ii], gamma_n[ii], T, U[,ii])
    })
    return(loglike)

}


#' Get loglike for SIR with U statistics
#' 
#' @param par c("beta", "gamma")
#' @param T number of steps 0:(T-1) inclusive
#' @param U  U a 3xN matrix with (A0, SMax, IMax) 
#' @param suff_stat_type character either "U" or "X" default is U
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


#' Get loglike for SIR with X statistics
#' 
#' @param par c("beta", "gamma")
#' @param T number of steps 0:(T-1) inclusive
#' @param suff_stat  a Tx3 matrix with (St, It, Rt).  Must have data at every point!!
#' @param suff_stat_type character either "U" or "A" default is U
#' @param fxn function either "KM" for Kermack and McKendrick prob function or "RF" for Reed-Frost prob
#' @param neg_loglike logical.  Default is TRUE.  Should we return the negative loglike?
#' @param use_exp_X logical.  Default is FALSE.  Should we use expected values of X for the prob function or the data?
#' @param x0 c(S0, I0, R0) or NULL
#' @param inf_nbrs TBD
#' @return scalar loglike
loglike_sir.X <- function(par, T, suff_stat, prob_fxn,
                          use_exp_X, x0,
                          inf_nbrs){

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
        X <- suff_stat
    }
    N <- sum(X[1,])
    pt <- fxn(par, X, inf_nbrs, N)

    loglike <- loglike_SIR_CM(T, pt, par[2], suff_stat, N)

    return(loglike)

}

#' Split U suff stat by time
#'
#' @param U a 3xN matrix with (A0, SMax, IMax)
#' @param tmin minimum time of last day prior to infection (inclusive).  Default is 0
#' @param tmax maximum time of last day prior to infection (inclusive).
#' @return subset of U such that t_min <= SMax <= t_max
split_U <- function(U, tmin = 0, tmax){

    ## either initially susceptible and last day infected between tmin and tmax OR
    ## initially infected and last time infected is minimum tmin
    inds <- which((U[1,] == 0 & tmin <= U[2, ] & U[2, ] <= tmax) |
                  U[1, ] == 1 & tmin == 0)
    new_U <- U[, inds]
    new_U[1, ] <- ifelse(new_U[2,] <= tmin, 1, new_U[1,]) # if last day susceptible is 
    new_U[2, ] <- ifelse(new_U[2,] <= tmin, tmax - tmin, new_U[2,] - tmin)
    new_U[3, ] <- ifelse(new_U[3,] <= tmin, tmax - tmin, new_U[3,] - tmin)
    return(new_U)
}



#' Update U after subsetting by time
#'
#' @param A0 initial states at time 0
#' @param SMax maximum time before infection (between tmin and tmax)
#' @param tmin
udpate_A0 <- function(A0, Smax, tmin, tmax){
    

}


#' Combine output from am_sir with different time ranges
#'
#' @param sims_X_list list of multiple "X" outputs from am_sir()
#' @param N total number of agents
#' @return a data frame where we add the corresponding states together
combine_sims_X <- function(sims_X_list, N){
    sims_X_list <- lapply(sims_X_list, function(df){
        if(class(df) != "data.frame"){
            df <- as.data.frame(df)
        }
        df
    })
    unique_times <- sort(unique(do.call('c', lapply(sims_X_list, function(df) df$t))))
    unique_sims <- sort(unique(do.call('c', lapply(sims_X_list, function(df) df$ll))))
    grid <- expand.grid(t = unique_times, ll = unique_sims)
    joined_list <- lapply(sims_X_list, function(df){
        new_X <- plyr::join(grid, df, by = c("t", "ll"))
        new_X$S <- ifelse(is.na(new_X$S), 0, new_X$S)
        new_X$I <- ifelse(is.na(new_X$I), 0, new_X$I)
        new_X$R <- ifelse(is.na(new_X$R), 0, new_X$R)
        ## Now adjust the later T
        T <- max(df$t, na.rm = TRUE)
        last_R <- df$R[T]
        last_I <- df$I[T]
        last_S <- df$S[T]
        new_X$R <- ifelse(new_X$t > T, last_R, new_X$R)
        new_X$I <- ifelse(new_X$t > T, last_I, new_X$I)
        new_X$S <- ifelse(new_X$t > T, last_S, new_X$S)
        return(new_X[, c("S", "I", "R")])
    })
    X <- Reduce('+', joined_list)
    df <- cbind(grid, X)
    df$S <- N - df$I - df$R
    return(df)
    
}
                      


#' Combine estimate of data frames perhaps at different times
#'
#' @param df_list list of estimate data frames
#' @param N total number of agents
#' @param CI whether we should include variance
#' @return combined data frame for all unique time steps
combine_ests <- function(df_list, N, CI = FALSE){
    L <- length(df_list)
    new_df_list <- vector(mode = "list", length = L)
    data_type <- df_list[[1]]$data_type[1]
    unique_times <- sort(unique(do.call('c', lapply(df_list, function(df) df$t))))
    out_df <- data.frame(t = unique_times, S_mean = 0, I_mean= 0, R_mean= 0)
    for(ii in 1:L){
        new_df <- data.frame(t = unique_times)
        new_df <- plyr::join(new_df, df_list[[ii]][, c("t", "S_mean", "I_mean", "R_mean")], by = "t")
        Tmax <- max(df_list[[ii]]$t, na.rm = TRUE)
        Tmin <- min(df_list[[ii]]$t, na.rm = TRUE)
        if(Tmax < max(unique_times)){
            new_df$S_mean <- ifelse(is.na(new_df$S_mean), new_df$S_mean[Tmax + 1], new_df$S_mean)
            new_df$I_mean <- ifelse(is.na(new_df$I_mean), new_df$I_mean[Tmax + 1], new_df$I_mean)
            new_df$R_mean <- ifelse(is.na(new_df$R_mean), new_df$R_mean[Tmax + 1], new_df$R_mean)
        }
        if(Tmin > min(unique_times)){
            new_df$S_mean <- ifelse(is.na(new_df$S_mean), new_df$S_mean[Tmin + 1] +
                                                          new_df$I_mean[Tmin + 1] +
                                                          new_df$R_mean[Tmin + 1], new_df$S_mean)
            new_df$I_mean <- ifelse(is.na(new_df$I_mean), 0, new_df$I_mean)
            new_df$R_mean <- ifelse(is.na(new_df$R_mean), 0, new_df$R_mean)
        }


        out_df$S_mean <- out_df$S_mean + new_df$S_mean
        out_df$I_mean <- out_df$I_mean + new_df$I_mean
        out_df$R_mean <- out_df$R_mean + new_df$R_mean
        ## max_R <- max(new_df$R_mean, na.rm = TRUE)
        ## new_df$R_mean <- ifelse(new_df$t > T, max_R, new_df$R_mean)
        ## out_df$I_mean <- new_df$I_mean + out_df$I_mean
        ## out_df$R_mean <- new_df$R_mean + out_df$R_mean
        
    }
##    out_df$S_mean <- N - out_df$I_mean - out_df$R_mean
    out_df$data_type <- data_type
    return(out_df)

}


#' Turn A0 into X matrix with G groups
#'
#' @param A0 vector of size N with initial states (0 is S, 1 is I, 2 is R)
#' @param G_id vector of length N where the entry is the group ID, groups should be 1 to G
#' @param G total number of groups
#' @param total X
#' @return a x0 vector where the first G entries are #S for each group, next G are I, and final G are R
sir_init_groups <- function(A0, G_id, G, X){
    if(is.null(G_id)){
        return(X[1,])
    }
    x0 <- numeric(G * 3)
    for(gg in unique(G_id)){
        x0[gg] <-  length(which(A0 == 0 & G_id == gg))
        x0[G + gg] <-  length(which(A0 == 1 & G_id == gg))
        x0[2 * G + gg] <-  length(which(A0 == 2 & G_id == gg))
    }
    return(x0)


}



#' Turn a stratified X_group into X
#'
#' @param X_group matrix of size T x (3*G) where G is the number of groups, S1, \dots, SG, I1, \dots IG, R1, \dots, RG is the order
#' @return Tx3 matrix where S = \sum_{g=1}^G Sg, and same for others
aggregate_X <- function(X_group){
    G <- ncol(X_group) / 3
    S <- rowSums(X_group[, 1:G, drop = FALSE])
    I <- rowSums(X_group[, (G+1):(2*G), drop = FALSE])
    R <- rowSums(X_group[, (2*G+1):(3*G), drop = FALSE])
    X <- cbind(S, I, R)
    colnames(X) <- NULL
    return(X)
}


#' Turn vector of size G into a vector of size N based on group assignments
#'
#' @param par vector of size G
#' @param G_id vector of size N with entries 1 to G
#' @param N total number of agents
#' @return par_n vector of size N
group_to_agent_par <- function(par, G_id, N){
    if(is.null(G_id)){
        par_n <- rep(par, N)
    } else{
        par_n <- par[G_id]
    }
    return(par_n)
}
    
