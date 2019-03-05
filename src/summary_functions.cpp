#include <Rcpp.h>
using namespace Rcpp;


// take 3xN suff stats and convert it to (Tmax)xN for SIR
// [[Rcpp::export]]
IntegerMatrix UtoA_SIR(IntegerMatrix U, int Tmax){
  int N = U.ncol();
  IntegerMatrix A(Tmax, N);
  // Fill in initial states
  for(int nn=0; nn < N; nn++){
    A(0, nn) = U(0, nn);
  } 
  // Initialize vars
  int SMaxTime = Tmax;
  int IMaxTime = Tmax;
  int A0 = 0;
  // Populate entries
  for(int nn=0; nn < N; nn++){
    A0 = U(0, nn);
    SMaxTime = U(1, nn);
    IMaxTime = U(2, nn);
    Rprintf("A0: %d, Smax: %d, IMax: %d\n", A0, SMaxTime, IMaxTime);
    for(int tt=1; tt < Tmax; tt++){
      if(A0 == 0){ // If agent starts susceptible
        if(SMaxTime < IMaxTime & SMaxTime < Tmax){  // agent gets infected and recovers
          if(tt <= SMaxTime){
            A(tt,nn) = 0;
          } else if (tt <= IMaxTime){
            A(tt,nn) = 1;
          } else{
            A(tt,nn) = 2;
          }
        } else if(SMaxTime < Tmax){ // agent gets infected but doesn't recover
          if(tt <= SMaxTime){
            A(tt,nn) = 0;
          } else{
            A(tt,nn) = 1;
          }
        } else {
          A(tt,nn) = 0; // agent doesn't get infected
        }
      } else if(A0 == 1){  // agent starts infected
        if(IMaxTime < Tmax){ // agent recovers
          if(tt <= IMaxTime){
            A(tt,nn) = 1;
          } else{
            A(tt,nn) = 2;
          }
        } else{ // agent does not recover
          A(tt, nn) = 1; // the rest are infectious
        }
      } else if(A0 == 2){ // agent starts recovered
        A(tt, nn) = 2; // all recovered
      }
    }
  }
  return A;
  
}

// A to U SIR
// [[Rcpp::export]]
IntegerMatrix AtoU_SIR(IntegerMatrix A){
  int Tmax = A.nrow();
  int N = A.ncol();
  IntegerMatrix U(3, N);
  // Initialize vars
  int A0;
  int SMax=Tmax;
  int IMax=Tmax;
  int Atn = -1;
  for(int nn=0; nn < N; nn++){
    SMax = Tmax -1;
    IMax = Tmax - 1;
    U(0,nn) = A(0, nn); // Initial state
    for(int tt=0; tt < Tmax; tt++){
     Atn = A(tt, nn); // current state
      if(Atn == 0){
        SMax = tt;
      } else if(Atn == 1){
        IMax = tt;
      }
    }
    U(1, nn) = SMax;
    U(2, nn) = IMax;
  }
  return U;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
U <- matrix(c(0, 0, 1, 0,
              0, 1, 2, 2,
              1, 2, 2, 2), byrow = TRUE, nrow = 3)
Tmax <- 3
A <- UtoA_SIR(U, Tmax)
A

A <- matrix(c(0, 0, 1, 2,
              0, 1, 1, 2,
              1, 2, 1, 2), byrow = TRUE, nrow = 3)
U <- AtoU_SIR(A)
U

UtoA_SIR(AtoU_SIR(A), Tmax)

U <- matrix(c(0, 0, 1, 0,
              0, 1, 2, 2,
              1, 2, 2, 2), byrow = TRUE, nrow = 3)
AtoU_SIR(UtoA_SIR(U, Tmax))

*/
