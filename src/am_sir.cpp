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

