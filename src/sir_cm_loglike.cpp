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




// compute the SIR loglike for the CM (where X is the sufficient statistic)
//
//@param T total number of steps (0:(T-1) inclusive
//@param pt vector of size T-1, entry i is the probability of becoming S from step i-1 to i for i=1, \dots, T-1
//@param gamma prob of recovery
//@param X Tx3 Matrix S, I, R
//@param N total number of agents
//@return loglike
// [[Rcpp::export]]
double loglike_SIR_CM(int T, NumericVector pt, double gamma, NumericMatrix X, int N) {

  double loglike_S = 0;
  double loglike_R = 0;
  double loglike = 0;
  double delta_S;
  double delta_R;
  for(int tt=1; tt < T; tt++){

    delta_S = 1.0 * X(tt-1, 0) - 1.0 * X(tt, 0);
  // Rprintf("delta_S %.1f\n", delta_S);
    delta_R = X(tt, 2) - X(tt-1, 2);
   // Rprintf("delta_R %.1f\n", delta_R);

    loglike_S = 0;
    loglike_R = 0;
    if( X(tt-1, 0) > 0){ // If old number of S is greater than 0
      loglike_S = delta_S * log(pt[tt-1]) + X(tt,0) * log(1 - pt[tt-1]);
    //  Rprintf("loglike_S %.5f\n", loglike_S);
     
    }
    if(X(tt-1, 1) > 0){ // If old number of I is greater than 0
      loglike_R = delta_R * log(gamma) + (X(tt-1,1) - delta_R) * log(1-gamma);
    //  Rprintf("loglike_R %.5f\n", loglike_R);
    }
    
    
    loglike = loglike + loglike_S + loglike_R;
    
  }
  return loglike;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
par <- c(.2, .7)
T <- 3
pt <- rep(par[1], T-1)
X <- matrix(c(3, 1, 0,
              2, 2, 0,
              1, 2, 1 ), byrow = TRUE, nrow = 3)
x <- X[1,]

N <- sum(x)


out <- loglike_SIR_CM(T, pt, par[2], X, N)
exp_loglike <- (1 * log(pt[1]) + 2 * log(1-pt[1])) +
  (0 + 1 * log(1-par[2])) +
  (1 * log(pt[2]) + 1 * log(1-pt[2])) +
  (1 * log(par[2]) + 1 * log(1-par[2]))
*/
