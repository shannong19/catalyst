#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int updateAgent_SIR(int current_state, double prob){
  int new_state = current_state;
  if(current_state == 2){
    return current_state;
  } else{
    new_state = current_state + as<int>(rbinom(1, 1, prob));
  }
  
  return new_state;
  
}


// [[Rcpp::export]]
IntegerVector AMTime_SIR(IntegerVector A,
                         NumericVector prob_inf, NumericVector prob_rec){
  int N = A.size();
  int current_state;
  double prob;
  IntegerVector newA(N);

  for(int nn = 0; nn < N; nn++){
    current_state = A[nn];
    if(current_state == 0){
      prob = prob_inf[nn];
    } else if(current_state == 1){
      prob = prob_rec[nn];
    } else {
      prob = 1.0;
    }

    newA[nn] = updateAgent_SIR(current_state, prob);
  }
  
  return newA;
}

// Get the number of neighbors who are infectious for each agent
// // Starts at 0!
// [[Rcpp::export]]
IntegerVector get_n_nbr_inf(List nbr_list, IntegerVector inf_inds){
  int N = nbr_list.size();
  IntegerVector nbr_inds;
  IntegerVector inf_nbr_inds;
  IntegerVector n_inf_nbrs(N);
  for(int nn=0; nn < N; nn++){
    nbr_inds = nbr_list[nn];
    inf_nbr_inds = intersect(inf_inds, nbr_inds);
  //  Rprintf("n inf nbrs: %d\n", inf_nbr_inds.size());
    n_inf_nbrs[nn] = inf_nbr_inds.size();
  }
  return n_inf_nbrs;
  
}



// Update the neighbor list based ongoing preventions
//
// @param prevention_types one of "isolation", "quarantine", or "school_closure"
// @param prevention_nbrs pre-computed list of nbrs post-prevention
// @param orig_nbrs original neighbor list
// @param inf_inds current infectious agent indices
// @param env_mat N x E data frame
// @param preventions_inds which agents will change their neighbors
// @details IMPORTANT:  INDEXING STARTS AT 0
// @return updated neighborl list
// [[Rcpp::export]]
List updateNeighbors(List orig_nbrs, 
                     List preventions_nbrs, 
                     IntegerVector preventions_inds){
  // Make new independent copy of original neighbors TODO: make more efficient
  List new_nbrs = clone(orig_nbrs);
  int n_inf = preventions_inds.size();
  if(n_inf == 0){
    return new_nbrs;
  }
  int index=0;
  for(int ii=0; ii < n_inf; ii++){
    index = preventions_inds[ii];
    new_nbrs[index] = preventions_nbrs[index];
  }
  
  return new_nbrs;
}


// Remove isolated/quarantined agents from the neighbors list
// 
// @param orign nbrs list of agents and their respective neighbors
// @param preventions_nbrs the indices of neighbors (starting at 0) who will have preventions implemented
// @param inf_inds indices of currently infectious individuals
// @param preventions_inds neighbors that are being removed from contact
// @return new_nbrs which has removed the prevention indices
// [[Rcpp::export]]
List removeNeighbors(List orig_nbrs, List preventions_nbrs,
                     IntegerVector inf_inds,
                     IntegerVector preventions_inds){
  // Make new independent copy of original neighbors TODO: make more efficient
  List new_nbrs = clone(orig_nbrs);
  IntegerVector inf_prev_inds;
  inf_prev_inds = intersect(inf_inds, preventions_inds);  // PRESUMES THIS IS SORTED
  int n_prev = inf_prev_inds.size();
  int n_set_diff = 0;
  if(n_prev == 0 |
     preventions_inds(0) == -1 | 
     inf_inds(0) == -1 
     ){  // no preventions to be made, return original neighbor list without change
    return new_nbrs;
  }
  // Loop through infected and preventions indices
  // Find their neighbors
  // Find their permanent neighbors
  // Take set difference of neighbors / perm neighbors
  // Loop over set differences indices
  // Remove the infected and prevention neighbor
  int prev_ind = 0;
  int nbr_ind = -1;
  IntegerVector all_nbr_inds;
  IntegerVector perm_nbr_inds;
  IntegerVector set_diff_inds;
  IntegerVector old_nbr_inds;
  for(int ii=0; ii < n_prev; ii++){
    prev_ind = inf_prev_inds(ii);
    all_nbr_inds = new_nbrs[prev_ind];
    //Rprintf("prev_ind: %d; all_nbr_inds(0): %d\n", prev_ind, all_nbr_inds(0));
    if(all_nbr_inds(0) == -1){ // if no neighbors, cannot infect skip to next
          continue;
    }
    perm_nbr_inds = preventions_nbrs[prev_ind];
    set_diff_inds = setdiff(all_nbr_inds, perm_nbr_inds);
    n_set_diff = set_diff_inds.size();
   // Rprintf("n_set_diff: %d\n", n_set_diff);
    if(n_set_diff == 0){ // if no neighbors to change, go to next
      continue;
    }
    for(int jj=0; jj < n_set_diff; jj++){  // get rid of the prev ind in the neighbor
      nbr_ind = set_diff_inds(jj);
      old_nbr_inds = new_nbrs[nbr_ind];
      if(old_nbr_inds.size() == 1){
        new_nbrs[nbr_ind] = preventions_nbrs[nbr_ind];
      } else {
        new_nbrs[nbr_ind] = old_nbr_inds[old_nbr_inds != prev_ind];
      }
    }
    new_nbrs[prev_ind] = perm_nbr_inds; // get rid of all nbrs for ind but the permanent ones 
  }
  
  return new_nbrs;

}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

nbr_list <- list(c(2, 3, 4), c(1, 3), c(1, 2), c(1), c(-1))
inf_inds <- c(1:4)

z <- get_n_nbr_inf(nbr_list, inf_inds)
z



current_state <- 1
prob <- 1
new_state <- updateAgent_SIR(current_state, prob)
new_state

## AMTime SIR
N <- 4
T <- 3
A <- matrix(c(0, 0, 0, 1,
              1, 0, 0, 1,
              2, 0, 1, 1), byrow = TRUE, nrow = T, ncol = N)
t <- 0
prob_inf <- rep(1, N)
prob_rec <- rep(1, N)
out <- AMTime_SIR(A[1,], prob_inf, prob_rec)
out
*/

