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



// Calculate the SIR  diff eqs for G groups
// 
// @param x initial values, first G are initial S, second G are initial I, and third G are initial R
// @param beta_vec vector of length G, values between 0 and 1
// @param gamma_vec vector of length G values between 0 and 1
// @param T integer 0:(T-1) inclusive (so length T)
// @param prob_type  0 is KM and 1 is RF
// @return matrix of Tx(3G)
// [[Rcpp::export]]
NumericMatrix sirLoopGroups(NumericVector x, NumericVector beta_vec, 
                      NumericVector gamma_vec, int T,
                      int prob_type){
  int K = x.size();
  int G = K / 3;
  double N = 0;
  NumericMatrix new_x(T, K);
 // Get total number in population and initialize new matrix
  for(int ii=0; ii < K; ii++){
    N += x[ii];
    new_x(0, ii) = x[ii];
  }
  // Get total number in each subpop
  NumericVector Nk(G);
  for(int ii=0; ii < G; ii++){
    Nk(ii) = x[ii] + x[ii + G] + x[ii + 2 * G];
  }



  for(int tt = 1; tt < T; tt++){
    double totalI = 0;  
    for(int ii=G; ii < G + G; ii++){
      totalI += new_x(tt-1, ii); // sum up the number of infectious individuals at previous time step
    }
    // Get new numbers of individuals in each group, where everyone has equal chance of interacting with ANY infected individual
    for(int gg=0; gg < G; gg++){
      if(prob_type == 1){  // KM
        new_x(tt, gg) = new_x(tt-1, gg) - beta_vec[gg] * new_x(tt-1, gg) * totalI / N;
      } else{ // RF
        new_x(tt, gg) = new_x(tt-1, gg) - 
          new_x(tt-1, gg)  * (1 - pow((1-beta_vec[gg]/N), totalI));
      }
      new_x(tt, 2 * G + gg) = new_x(tt-1, 2 * G + gg) + gamma_vec[gg] * new_x(tt-1, G + gg);
      new_x(tt, G + gg) = Nk[gg] - new_x(tt, gg) - new_x(tt, 2 * G + gg);
      
    
    }
  }
  
  return new_x;
}






// for prob_type 0 is KM and 1 is RF
// [[Rcpp::export]]
NumericVector sirLoop(NumericVector x, double beta, double gamma, int T,
                      int prob_type){
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

// Make list of neighbors
// 
// @param env_mat matrix of size N x E0 if has no assignment
// @return list where each entry is a column of integers.  Each integer is a neighbor.  If it is -1, then the agent has no neighbors.
// VERY IMPORTANT:  INDEXING FOR NEIGHBORS STARTS AT 0
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
           // Rprintf("n: %d, m %d, ee %d \n", nn, mm, ee);
            nbr_vec[ii] = mm;
            ii = ii + 1;
            //Rprintf("nbr vec: %d \n", nbr_vec[ii]);
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

outG <- sirLoopGroups(x, beta, gamma, T, 0)
tail(outG)


out <- sirLoop(x, beta, gamma, T, 0)
tail(out)
head(out)
testthat::expect_equal(outG, out)
## Try multiple groups
x <- c(400, 550, 50, 0, 0, 0)
beta_vec <- c(.3, .1)
gamma_vec <- c(.1, .1)
T <- 5
outG <- sirLoopGroups(x, beta_vec, gamma_vec, T, 0)
outG
testthat::expect_true(all(rowSums(outG[, c(2, 4, 6)]) == 550))


env_mat <- matrix(c(0, 1, 1,
                    1, 2, 1,
                    0, 1, 1,
                    0, 0, 2), byrow = TRUE, nrow = 4)

map <- initializeNeighbors(env_mat)
map
***/


