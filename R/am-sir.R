## SKG
## March 5, 2019
## V5? of catalyst AM.



#' Run simulations of AM SIR
#' 
#' @param LL total number of simulations to run.  Default is 1
#' @param T total amount of time steps (0:(T-1))
#' @param A0 vector of length N of intial values for N agents.  0- S, 1 - I, 2 - R
#' @param prob_fxn character either "KM" for Kermack and McKendrick diff eqs. or "RF" for Reed Frost chain binomial. See details for more info
#' @param par1_vec vector of length N ranging between 0 and 1.  This corresponds to "beta" for KM and "rho" for RF
#' @param gamma_vec vector of length N ranging between 0 and 1.  Corresonds to "gamma" the average inverse time to recovery
#' @param nbr_list list/map where the key is the is the index of the agent and values are its neigbhors or -1  if it has no neighbors.  If nbr_list = NULL, we assume everyone is a neighbor to each other
#' @param use_exp_X default is TRUE.  If not TRUE, it will use previous values from simulation
#' @param keep_A logical.  Default is FALSE
#' @param keep_U logical.  Default is FALSE
#' @param write_sim logical. Default is FALSE
#' @param writing_list list of writing parameters for write_sim.  Default is NULL.
#' @param do_par logical.  Default is FALSE.  Coming soon
#' @param do_preventions logical.  Default is FALSE
#' @param preventions_type one of "isolation" or "quarantine".  Default is NULL
#' @param preventions_nbrs neighbor list of pre-comptued neighbors given the prevention is in place.  Default is NULL
#' @param env_df NxE data frame.  Default is NULL
#' @param preventions_vars variable names which we subset to when applying preventions.  eg isolation means only neighbors are household ID.
#' @return list of X a Tx3 matrix of number of S, I, and R at each time step;  A a TxN matrix of each agents state at each path or NULL if keep_A is FALSE
am_sir <- function(L, T,
                   A0, prob_fxn = "KM",
                   par1_vec, gamma_vec,
                   nbr_list = NULL,
                   use_exp_X = TRUE,
                   keep_A = FALSE,
                   keep_U = FALSE,
                   write_sim = FALSE,
                   writing_list = NULL,
                   do_par = FALSE,
                   do_preventions = FALSE,
                   preventions_type = NULL,
                   preventions_nbrs = NULL,
                   env_mat = NULL,
                   preventions_vars = NULL){

    ## Initialize
    N <- length(A0)
    A <- matrix(0, nrow = T, ncol = N)
    A[1,] <- A0
    

    if(use_exp_X){
      pf <- ifelse(prob_fxn == "KM", 0, 1)
      X0 <- AtoX(A[1,, drop = FALSE], 3)
      X_theor <- sirLoop(X0, par1_vec[1], gamma_vec[1], T, pf)
    } else{
        X_theor <- NULL
    }

    if(!do_par){
        X_list <- vector(mode = "list", length = L)
        
        if(keep_A){
            A_arr <- array(0, dim = c(L, T, N))
        } else{
            A_arr = NULL
        }
        if(keep_U){
            U_arr <- array(0, dim = c(L, 3, N))
        } else{
            U_arr <- NULL
        }

        for(ll in 1:L){
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
                      writing_list = writing_list,
                      do_preventions = do_preventions,
                      preventions_type = preventions_type,
                      prevention_nbrs = preventions_nbrs,
                      env_df = env_df,
                      preventions_vars = preventions_vars)
              X <- out$X
              X <- cbind(c(0:(T-1)), X, rep(ll, T))
              colnames(X) <- c("t", "S", "I", "R", "ll")
              X_list[[ll]] <- X

              ## U and A
              if(keep_A){
                  A_arr[ll,,] <- out$A
              }
              if(keep_U){
                  U_arr[ll,,] <- out$U
              }
              
        }

    } else{
        stop("no || option at the moment")
    }
    
    return(list(X_mat = do.call('rbind', X_list), A = A_arr, U = U_arr))

}


#' Run one simulation of AM SIR
#' 
#' @param ll simulation number.  Default is 1
#' @param A TxN matrix where first row has correct intial values.  0- S, 1 - I, 2 - R
#' @param prob_fxn character either "KM" for Kermack and McKendrick diff eqs. or "RF" for Reed Frost chain binomial. See details for more info
#' @param par1_vec vector of length N ranging between 0 and 1.  This corresponds to "beta" for KM and "rho" for RF
#' @param gamma_vec vector of length N ranging between 0 and 1.  Corresonds to "gamma" the average inverse time to recovery
#' @param nbr_list list/map where the key is the is the index of the agent and values are its neigbhors or -1  if it has no neighbors.  If nbr_list = NULL, we assume everyone is a neighbor to each other
#' @param X_theor default is NULL.  If not null, it is the expected/theoretical values of X at each time step
#' @param keep_A logical.  Default is FALSE
#' @param keep_U logical.  Default is FALSE
#' @param write_sim logical. Default is FALSE
#' @param writing_list list of writing parameters for write_sim.  Default is NULL.
#' @param do_preventions logical.  Default is FALSE
#' @param preventions_type one of "isolation" or "quarantine".  Default is NULL
#' @param prevention_nbrs neighbor list of pre-comptued neighbors given the prevention is in place.  Default is NUL
#' @param env_df NxE data frame.  Default is NULL
#' @param preventions_vars variable names which we subset to when applying preventions.  eg isolation means only neighbors are household ID.
#' @return list of X a Tx3 matrix of number of S, I, and R at each time step;  A a TxN matrix of each agents state at each path or NULL if keep_A is FALSE
am_sir_one_sim <- function(ll = 1, A, prob_fxn = "KM",
                           par1_vec, gamma_vec,
                           nbr_list = NULL,
                           X_theor = NULL,
                           keep_A = FALSE,
                           keep_U = FALSE,
                           write_sim = FALSE,
                           writing_list = NULL,
                           do_preventions = FALSE,
                           preventions_type = NULL,
                           env_df = NULL,
                           preventions_vars = NULL){

    U <- NULL
    T <- nrow(A)
    N <- ncol(A)

    out_A <- NULL
    orig_nbrs <- nbr_list
    if(do_preventions){  ## pre-compute prevention_neighbors
        if(is.null(nbr_list)){
            ## Populate neighbors list
            mat <- matrix(rep(1, N), ncol = 1)
            nbr_list <- initializeNeighbors(mat)
        }
        preventions_nbrs <- precompute_preventions_nbrs(env_df, preventions_vars)
    }
    for(tt in 1:(T-1)){
        if(do_preventions){
            nbr_list <-updateNeighbors(orig_nbrs, preventions_nbrs, inf_inds)
        }
        inf_inds <- which(A[tt,] == 1)
        prob_inf <- get_prob_inf(par1_vec, inf_inds,
                                 nbr_list,
                                 X_theor[tt,],
                                 prob_fxn)
        A[tt+1, ] <- AMTime_SIR(A[tt,], prob_inf, gamma_vec)

    }
    X <- AtoX(A, 3)
    if(keep_U){
        U <- AtoU_SIR(A)
    }
    if(!keep_A){
        A <- NULL
    }
    
    
    return(list(X = X, U = U, A = A))



}

#' @param par1_vec vector of length N ranging between 0 and 1.  This corresponds to "beta" for KM and "rho" for RF
#' @param inf_inds indices of infectious neigbhors
#' @param nbr_list list/map where the key is the is the index of the agent and values are its neigbhors or -1  if it has no neighbors.  If nbr_list = NULL, we assume everyone is a neighbor to each other
#' @param X_theor default is NULL.  If not null, it is the expected/theoretical values of X at each time step 
#' @param prob_fxn character either "KM" for Kermack and McKendrick diff eqs. or "RF" for Reed Frost chain binomial. See details for more info
#' @return prob_inf vector of length N, which the probability of current agent getting infected
get_prob_inf <- function(par1_vec, inf_inds,
                         nbr_list = NULL,
                         X_theor = NULL,
                         prob_fxn = "KM"){

    N <- length(par1_vec)
    if(!is.null(X_theor)){ # Use expected X values
        n_inf_nbrs <- X_theor[2]
    } else if(is.null(nbr_list)){
        nbr_inf_inds <- inf_inds
        n_inf_nbrs <- length(nbr_inf_inds)
    } else {
            n_inf_nbrs <- get_n_nbr_inf(nbr_list, inf_inds)
    }
    
    if(prob_fxn == "KM"){
        prob_inf <- par1_vec * n_inf_nbrs / N
    } else if(prob_fxn == "RF"){
         prob_inf <- 1 - (1 - par1_vec / N)^n_inf_nbrs
    } else{
        stop("That probability function doesn't exist")
    }
    return(prob_inf)
}
