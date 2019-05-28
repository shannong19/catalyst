#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
int findIfSus(int A0, int IMax, int T){
  // agent is susceptible only IF:
  // A0 is 0
  // AND
  // (IMax = T-1 )
  int isSus = 0; // 1 is susceptible
  if(A0 == 0 & IMax == T-1){
    isSus = 1;
  }
  return isSus;
  
}


// [[Rcpp::export]]
IntegerVector whichState(IntegerVector x, int state){
  int N = x.size();
  IntegerVector inds(N, -1);
  int counter = 0;
  for(int ii=0; ii < N; ii++){
    if(x[ii] == state){
      inds[ii] = ii; 
      counter++;
    }
  }
  if(counter == 0){ // nothing to return
    IntegerVector inds(1);
    inds[0] = -1;
    return inds;
  }
  return inds[inds > -1];
}

// Quick and dirty AM SIR
// 
// Return the U matrix for a single iteration
// [[Rcpp::export]]
IntegerMatrix AMSIR_inf_inner(int ll, int T,
                    IntegerVector A0,
                    List nbrList,
                    double beta,
                    double gamma){
  
  int N = A0.size();
  IntegerMatrix U(3, N);
  IntegerVector SMax(N, T-1);  // Default is T-1 the max value
  IntegerVector IMax(N, T-1); // Same for max I
  IntegerVector infInds = whichState(A0, 1);
  IntegerVector newA0 = clone(A0);
  int infInd;
  int isSus;
  int nbrInd;
  int isNewInf=0;
  int isNewRec=0;
  double infProb;
  for(int tt=1; tt < (T-1); tt++){
    infInds = whichState(newA0, 1);
    // Loop through infectious
    if(infInds[0] = -1){ // If there are no infectious, stop
      break;
    }
    for(int ii=0; ii < infInds.size(); ii++){
      infInd = infInds[ii];
      IntegerVector nbrVec;
      nbrVec = nbrList[infInd];
      if(nbrVec[0] = -1){
        continue;
      }
      // Loop through neighbors of infInd
      for(int jj=0; jj < nbrVec.size(); jj++){
        // check if susceptible
        nbrInd = nbrVec[jj];
        isSus = findIfSus(A0[nbrInd], IMax[nbrInd], T);
        if(isSus == 1){
          // Infect
          infProb = beta / (1.0 * nbrVec.size()); // Warning:  this is a big change! (from N)
          isNewInf = as<int>(rbinom(1, 1, infProb));
          if(isNewInf == 1){
            SMax[nbrInd] = tt-1;
            newA0[nbrInd] = 1; // tell model agent is infectious for next step
          }
         
        }
      }
      // Recover the infectious
      isNewRec = as<int>(rbinom(1, 1, gamma));
      if(isNewRec == 1){
        IMax[infInd] = tt-1;
        newA0[infInd] = 0; // tell model agent is no longer infectious
      }
      
    }
    
  }
  U(0, _ ) = A0;
  U(1, _ ) = SMax;
  U(2, _ ) = IMax;
  
  
  return U;
  

  
}

// LoopType is 0 for looping over infectious and 1 for over susceptible (but smartly, looking at who was infected)
// Return the list of U arrays for each iteration L
// [[Rcpp::export]]
IntegerMatrix AMSIR(int L, int T,
                              IntegerVector A0,
                              List nbrList,
                              double beta,
                              double gamma,
                              int loopType){
  
  List UList(L);
  for(int ll=0; ll < L; ll++){
    Ulist[ll] = AMSIR_inf_inner(ll, T,
                                A0, nbrList,
                                beta, gamma);
  }
  return Ulist;
  
} 




/*** R
x <- c(1, 0, 0)
out <- whichState(x, 0)
out

*/
