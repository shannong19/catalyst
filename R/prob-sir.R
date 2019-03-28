

#' Get the Kermack McKendrick probabilities of becoming infectious
#'
#' @param par beta, gamma or rho, gamma
#' @param X Tx3 matrix of S, I, R
#' @param inf_nbrs TBD.  Default is NULL
#' @param N total number of agents
#' @return vector of size T-1 where entry i is probability of getting infected from i-1 to i for i=1, \dots, T-1
KM_prob <- function(par, X, inf_nbrs = NULL, N){
    T <- nrow(X)
    pt <- par[1] * X[1:(T-1), 2] / N
    return(pt)
}


#' Get the Reed Frost probabilities of becoming infectious
#'
#' @param par beta, gamma or rho, gamma
#' @param X Tx3 matrix of S, I, R
#' @param inf_nbrs TBD.  Default is NULL
#' @param N total number of agents
#' @return vector of size T-1 where entry i is probability of getting infected from i-1 to i for i=1, \dots, T-1
RF_prob <- function(par, X, inf_nbrs = NULL, N){
    T <- nrow(X)
    pt <- 1 - (1 - par[1] / N)^X[1:(T-1), 2]
    return(pt)
}



#' Get the Kermack McKendrick probabilities of becoming infectious for groups
#'
#' @param params vector of betas of length N
#' @param X Tx3 matrix of S, I, R
#' @param inf_nbrs TBD.  Default is NULL
#' @param N total number of agents
#' @return  matrix of dimension  (T-1) X N where entry in is probability of an agent n of getting infected from i-1 to i for i=1, \dots, T-1
KM_probG <- function(params, X, inf_nbrs = NULL, N){
    T <- nrow(X)
    pt <- do.call('cbind', lapply(params, function(x){
        x * X[1:(T-1), 2] / N
        }))
    return(pt)
}


#' Get the Reed Frost probabilities of becoming infectious for groups
#'
#' @param params vector of betas of length N
#' @param X Tx3 matrix of S, I, R
#' @param inf_nbrs TBD.  Default is NULL
#' @param N total number of agents
#' @return  matrix of dimension  (T-1) X N where entry in is probability of an agent n of getting infected from i-1 to i for i=1, \dots, T-1
RF_probG <- function(params, X, inf_nbrs = NULL, N){
    T <- nrow(X)
    pt <- do.call('cbind', lapply(params, function(x){
        1 - (1 - x/ N)^X[1:(T-1), 2]
    }))
    return(pt)
}
