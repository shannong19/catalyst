#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


// for prob_type 0 is KM and 1 is RF
// [[Rcpp::export]]
NumericVector sirLoop(NumericVector x, double beta, double gamma, int T,
                      int prob_type) {
  double N = x[0] + x[1] + x[2];
  NumericMatrix new_x(T, 3);
  // Initialize matrix
  new_x(0, 0) = x[0];
  new_x(0, 1) = x[1];
  new_x(0, 2) = x[2];

  for(int tt = 1; tt < T; tt++){
    if(prob_type == 1){
      new_x(tt, 0) = new_x(tt-1, 0) - beta * new_x(tt-1,0) * new_x(tt-1, 1) / N;
      new_x(tt, 2) = new_x(tt-1, 2) + gamma * new_x(tt-1, 1);
      new_x(tt, 1) = N - new_x(tt, 0) - new_x(tt, 2);
    }
    else{
      new_x(tt, 0) = new_x(tt-1, 0) - 
        new_x(tt-1,0)  * (1 - pow((1-beta/N), new_x(tt-1, 1)));
      new_x(tt, 2) = new_x(tt-1, 2) + gamma * new_x(tt-1, 1);
      new_x(tt, 1) = N - new_x(tt, 0) - new_x(tt, 2);
    }
  }

  return new_x;
} 

// [[Rcpp::export]]
std::map<int, IntegerVector> initializeNeighbors(IntegerMatrix env_mat){
  std::map<int, IntegerVector> nbr_map;
  // iterator and pair
  std::map<int, IntegerVector>::iterator it;
  // pair <map<int, IntegerVector>::iterator, bool> ptr;

  int N = env_mat.nrow();
  int E = env_mat.ncol();
  IntegerVector nbr_vec(N);
  int ref_env = 0;
  int nbr_env = 0;
  for(int nn=0; nn < N; nn++){
    int ii = 0;
    for(int mm=0; mm < N; mm++){
       for(int ee=0; ee < E; ee++){
        ref_env = env_mat(nn, ee);
        if((ref_env != 0) & (nn != mm)){
          nbr_env = env_mat(mm, ee);
          if(ref_env == nbr_env){
            Rprintf("n: %d, m %d, ee %d \n", nn, mm, ee);
            nbr_vec[ii] = mm;
            ii = ii + 1;
            Rprintf("nbr vec: %d \n", nbr_vec[ii]);
            break;
          }
        }
       }

    }
    if(ii == 0){
      nbr_vec[ii] = -1;
      ii = 1;
    }
    // Shorten the vector
    IntegerVector assign_vec(ii);
    for(int jj=0; jj < ii; jj++){
      assign_vec[jj] = nbr_vec[jj];
    }

    // Add to map
    //nbr_map.insert(std::pair<int, IntegerVector>(nn, assign_vec));
    nbr_map[nn] = assign_vec;
  }


  return nbr_map;


}








// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
x <- c(950, 50, 0)
beta <- .1
gamma <- .03
N <- sum(x)
T <- 100

out <- sirLoop(x, beta, gamma, T, 0)
head(out)

env_mat <- matrix(c(0, 1, 1,
                    1, 2, 1,
                    0, 1, 1,
                    0, 0, 2), byrow = TRUE, nrow = 4)

map <- initializeNeighbors(env_mat)
map
***/


