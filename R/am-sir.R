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
#' @param preventions_type one of "isolation", "quarantine", or "closure".  Default is NULL
#' @param env_df NxE data frame.  Default is NULL
#' @param preventions_vars variable names which we subset to when applying preventions.  eg isolation means only neighbors are household ID.
#' @param days_infectious vector of size N.  number of days agent has been infected
#' @param preventions_delay  how many days we wait to do preventions
#' @param do_closure logical. Default is FALSE
#' @param closure_thresh list of length closure_vars where each entry of the list is the closure threshold percentage for the corresponding var
#' @param closure_vars which variables in env_df we will close on
#' @param closure_inds_list list of lists.  first level corresponds to the variable which we are looking at.  second level gives the category of hte variable we are looking at.  in the second level list, there are the indices of all the agents who belong to the same variable category
#' @param closure_time list of length closure_vars where each entry of the list is the time we will close the structure once closure is implemented
#' @param closure_max_T list of length closure_vars where each entry of the list is the max time we will close the structure
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
                   env_df = NULL,
                   preventions_vars = NULL,
                   days_infectious = NULL,
                   preventions_delay = 0,
                   do_closure = FALSE,
                   closure_thresh = list(.5),
                   closure_vars = NULL,
                   closure_inds_list = NULL,
                   closure_time = list(0),
                   closure_max_T = list(1)
                   ){


    print(paste("delay:", preventions_delay))
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

            if((ll %% 100) == 0){
                print(ll)
            }
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
                                  env_df = env_df,
                                  preventions_vars = preventions_vars,
                                  days_infectious = days_infectious,
                                  preventions_delay = preventions_delay,
                                  do_closure = do_closure,
                                  closure_thresh = closure_thresh,
                                  closure_vars = closure_vars,
                                  closure_inds_list = closure_inds_list,
                                  closure_time = closure_time,
                                  closure_max_T = closure_max_T)
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
#' @param env_df NxE data frame.  Default is NULL
#' @param preventions_vars variable names which we subset to when applying preventions.  eg isolation means only neighbors are household ID.  Default is NULL
#' @param days_infectious vector of size N.  number of days agent has been infected
#' @param preventions_delay  how many days we wait to do preventions
#' @param closure_thresh list of length closure_vars where each entry of the list is the closure threshold percentage for the corresponding var
#' @param closure_vars which variables in env_df we will close on
#' @param closure_inds_list list of lists.  first level corresponds to the variable which we are looking at.  second level gives the category of hte variable we are looking at.  in the second level list, there are the indices of all the agents who belong to the same variable category
#' @param closure_time list of length closure_vars where each entry of the list is the time we will close the structure once closure is implemented
#' @param closure_max_T list of length closure_vars where each entry of the list is the max time we will close the structure
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
                           preventions_vars = NULL,
                           days_infectious = NULL,
                           preventions_delay = 0,
                           do_closure = FALSE,
                           closure_thresh = NULL,
                           closure_vars = NULL,
                           closure_inds_list = NULL,
                           closure_time = list(0),
                           closure_max_T = list(1)){

    U <- NULL
    T <- nrow(A)
    N <- ncol(A)

    out_A <- NULL

##f    preventions_nbrs <- NULL
    if(do_preventions | do_closure){  ## pre-compute prevention_neighbors
        if(is.null(nbr_list)){
            ## Populate neighbors list
            mat <- matrix(rep(1, N), ncol = 1)
            nbr_list <- initializeNeighbors(mat)  ## This assumes EVERYONE is a neighbor at first
        }
        if(do_preventions){
            preventions_nbrs <- precompute_preventions_nbrs(env_df, preventions_vars) # these are the nbrs that are permanent
            preventions_inds <- 1:length(preventions_nbrs) -1
        }
    }

    
    orig_nbrs <- nbr_list

    for(tt in 1:(T-1)){
   ##     print(paste("tt", tt))
        inf_inds <- which(A[tt,] == 1) - 1 ## to start at index 0
        if(length(inf_inds) == 0 | inf_inds[1] == -1){  ## if there are no infectious left, there is no more movement, can output df
            for(ss in tt:(T -1)){
                A[ss + 1, ] <-  A[tt, ]
            }
            break
        }
        ##  print(paste("inf inds:", inf_inds))
        ## Isolation and quarantine
        if(do_preventions){
            browser()
            preventions_list <- update_preventions(nbr_list,
                               preventions_nbrs,
                               preventions_inds,
                               inf_inds,
                               days_infectious,
                               preventions_delay,
                               preventions_type)
            nbr_list <- preventions_list$new_nbrs
            preventions_inds <- preventions_list$new_preventions_inds
            do_preventions <- preventions_list$new_do_preventions
            days_infectious <- preventions_list$new_days_infectious
        }
        ## School closures
        if(do_closure){
            closure_list <- update_closure(nbr_list,
                                           inf_inds,
                                           closure_thresh,
                                           closure_vars,
                                           closure_inds_list,
                                           closure_time,
                                           closure_max_T)
            do_closure <- closure_list$do_closure
            nbr_list <- closure_list$new_nbrs
            closure_time <- closure_list$new_closure_time
            
        }
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


#' Update preventions
#'
#' @param nbr_list list/map where the key is the is the index of the agent and values are its neigbhors or -1  if it has no neighbors.  If nbr_list = NULL, we assume everyone is a neighbor to each other
update_preventions <- function(nbr_list,
                               preventions_nbrs,
                               preventions_inds,
                               inf_inds,
                               days_infectious = NULL,
                               preventions_delay = 0,
                               preventions_type = "isolation"){
    if(preventions_type %in% c("isolation", "quarantine")){
        preventions_list <- update_preventions.quar_iso(
            nbr_list,
            preventions_nbrs,
            preventions_inds,
            inf_inds,
            days_infectious,
            preventions_delay,
            preventions_type)
    }
    return(preventions_list)
}




#' Update preventions
#'
#' @param nbr_list
#' @param preventions_nbrs list of permanent neighbors
#' @param preventions_inds vector of indices we will do preventions on.  indexing starts at 0
#' @param inf_inds indices of infectious agents.  indexing starts at 0
#' @param days_infectious vector of size N.  number of days agent has been infected
#' @param preventions_delay  how many days we wait to do preventiosn
#' @return list with new nbr_list, prevention_inds, and whether we should continue preventions
update_preventions.quar_iso <- function(nbr_list,
                               preventions_nbrs,
                               preventions_inds,
                               inf_inds,
                               days_infectious = NULL,
                               preventions_delay = 0,
                               preventions_type = "isolation"){


   
    do_preventions <- TRUE
    ## If nothing to prevent, change do preventions and return the rest
    if(length(inf_inds) == 0 | inf_inds[1] == -1){
        do_preventions <- FALSE  ## no need to prevent if 0 infectious left


        return(list(new_nbrs = nbr_list,
                    new_preventions_inds = preventions_inds,
                    new_do_preventions = do_preventions,
                    new_days_infectious = days_infectious))
    }

    ## Otherwise check which infectious we apply preventions to
    if(!is.null(days_infectious)){
        if(length(inf_inds) > 0 & inf_inds[1] != -1){
            ## Add one to the days infectious
            if(preventions_type == "isolation"){
                days_infectious[inf_inds + 1] <- days_infectious[inf_inds + 1] + 1
            } else if(preventions_type == "quarantine"){
                quar_inds <- get_quar_inds(inf_inds, preventions_nbrs)
                days_infectious[quar_inds + 1] <- days_infectious[quar_inds + 1] + 1
            }
        }
        preventions_logical <- ifelse(days_infectious == preventions_delay,
                                      1, 0)  ## Only add NEW cases to the preventions list.
        preventions_inds <- which(preventions_logical == 1) - 1
        if(length(preventions_inds) == 0){
            preventions_inds <- -1
        }

    }


      ## remove infectious neighbors who have preventions implemented on them from the neighbors list.  Keep the preventions nbrs.  These people will always be susceptible.
    ## Remaining neighbors are unchanged.
    ##TODO:  inefficient because this makes new copy of nbr_list every time step  
    nbr_list <- removeNeighbors(nbr_list, preventions_nbrs,
                                preventions_inds, inf_inds)

    ## Update preventions_inds if days_infectious is null
    if(is.null(days_infectious)){
        preventions_inds <- preventions_inds[-which(preventions_inds %in%
                                                   inf_inds)]
        if(length(preventions_inds) == 0){
            preventions_inds <- -1
        }
    }

    return(list(new_nbrs = nbr_list,
                new_preventions_inds = preventions_inds,
                new_do_preventions = do_preventions,
                new_days_infectious = days_infectious))
}

#' Extract the indices to quarantine
#'
#' @param inf_inds currently infectious.  indexing starts at 0.
#' @param nbr_list list of neighbors
#' @return quarantine inds.  indexing starts at 0!
get_quar_inds <- function(inf_inds, nbr_list){
    nbr_inds <- do.call('c', nbr_list[inf_inds + 1])
    quar_inds <- sort(unique(c(nbr_inds, inf_inds)))
    return(quar_inds)
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
        ## Change proportion o infectious
       # n_nbrs <- get_n_nbr_inf(nbr_list, 1:N - 1)
       # N <- ifelse(n_nbrs == 0, 1, n_nbrs) ## avoid divide by 0
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
