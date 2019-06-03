#include <Rcpp.h>
using namespace Rcpp;



// // [[Rcpp::export]]
// int findIfSus(int A0n, int IMax, int T){
//   // agent is susceptible only IF:
//   // A0 is 0
//   // AND
//   // (IMax = T-1 )
//   int isSus = 0; // 1 is susceptible
//   if(A0n == 0 & IMax == T-1){
//     isSus = 1;
//   }
//   return isSus;
//   
// }


// [[Rcpp::export]]
IntegerVector makeStateVec(IntegerVector A0, int state){
  int N = A0.size();
  IntegerVector stateVec(N, -1);
  for(int ii=0; ii < N; ii++){
    if(A0[ii] == state){
      stateVec[ii] = ii;
    }
    
  }
  return stateVec;
  
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


// [[Rcpp::export]]
int countIntersect(IntegerVector infVec, IntegerVector nbrOfSusInds){
  int nNbr = nbrOfSusInds.size();
  //Rf_PrintValue(nbrOfSusInds);
  int count = 0;
  int ind = 0;
  for(int ii=0; ii < nNbr; ii++){
    ind = nbrOfSusInds[ii];
    if(infVec[ind] > -1){
      count++;
    }
  }
  return count;
}


// Take old vector keeping track of all the agent states 
// and make the indices matching inds and make them a positive value
// [[Rcpp::export]]
IntegerVector updateStateByInds(IntegerVector vec, IntegerVector inds){
  IntegerVector newVec = clone(vec);
  int nInds = inds.size();
  int ind=0;
  for(int ii=0; ii < nInds; ii++){
    ind = inds[ii];
    newVec[ind] = ind;
  }
  
  return newVec;
}

// [[Rcpp::export]]
IntegerVector updateStateVec(IntegerVector vec1, IntegerVector vec2){
  int N = vec1.size();
  IntegerVector newVec(N, -1);
  for(int ii=0; ii < N; ii++){
    if(vec1[ii] > -1 & vec2[ii] > -1){
      newVec[ii] = ii;
    }
  }
  return newVec;
}


// Quick and dirty AM SIR
// 
// Loop over first the susceptibles (subsetting to only neighbors of infectious)
// Return the U matrix for a single iteration
// [[Rcpp::export]]
IntegerMatrix AMSIR_sus_inner(int ll, int T,
                              IntegerVector A0,
                              List nbrList,
                              double beta,
                              double gamma){
  
  IntegerVector newA0 = clone(A0);
  int N = newA0.size();
  IntegerMatrix U(3, N);
  IntegerVector SMax(N, T-1);  // Default is T-1 the max value
  IntegerVector IMax(N, T-1); // Same for max I
  IntegerVector infInds;
  IntegerVector susInds;
  IntegerVector nbrSusInds;
  IntegerVector nbrInds;
  IntegerVector susVec(N, -1);
  IntegerVector infVec(N, -1);
  IntegerVector nbrSusVec(N, -1);
  IntegerVector nbrVec(N, -1);
  IntegerVector curNbrInds;
  IntegerVector nbrOfSusInds;
  int nInfNbrs;
  int infInd;
  int nbrSusInd;
  int nInfInds;
  int nNbrs = 1;
  int isNewRec;
  int isNewInf;
  double infProb;
  IntegerVector currentInfVec;
  //
  //
  susVec = makeStateVec(newA0, 0);
  infVec = makeStateVec(newA0, 1);
  for(int tt=1; tt < (T-1); tt++){
    IntegerVector nbrVec(N, -1); // reset to o neighbors
    infInds = infVec[infVec > -1];
    susInds = susVec[susVec > -1];
   // Rprintf("tt %d\n", tt);
   // Rprintf("size of infInds %d\n", infInds.size());
   // Rf_PrintValue(infInds);
    currentInfVec = clone(infVec);

    if(infInds.size() == 0){ // If there are no infectious, stop
     // Rprintf("no infectious inds at time t=%d", tt);
      break;
    }
    for(int ii=0; ii < infInds.size(); ii++){
      // Find susceptible neighbors, and combine them together
      infInd = infInds[ii];

      curNbrInds = nbrList[infInd];
      if(curNbrInds.size() > 0 & curNbrInds[0] != -1){
        nbrVec = updateStateByInds(nbrVec, curNbrInds); // this is the the unioned neighbors of all the infectious
      }

      // Recover the infectious
      isNewRec = as<int>(rbinom(1, 1, gamma));
      //   Rprintf("infInd: %d; isNewRec: %d; tt-1: %d\n", infInd, isNewRec, tt-1);
      if(isNewRec == 1){
        IMax[infInd] = tt-1;
        infVec[infInd] = -1; // tell model agent is no longer infectious
      }
   //   Rf_PrintValue(infInds);
      
    }
     nbrSusVec = updateStateVec(susVec, nbrVec);
     nbrSusInds = nbrSusVec[nbrSusVec > -1];
  //  Rprintf("size of nbrSusVec %d \n", nbrSusInds.size());
    if(nbrSusInds.size() != 0){
    //  Rprintf("drawing new infections");
      for(int jj=0; jj < nbrSusInds.size(); jj++){
        // Loop over the susceptible neighbors of the infectious
        nbrSusInd = nbrSusInds[jj]; // Should have at least one infectious neighbor
        nbrOfSusInds = nbrList[nbrSusInd];  // extract all their neighbors
        nNbrs = nbrOfSusInds.size();
        //Rprintf("nNbrs %d\n", nNbrs);
      // count total number
        if(nNbrs >= (N-1)){
          nInfInds = infInds.size();
        }else {
          nInfInds = countIntersect(currentInfVec, nbrOfSusInds); // count how many are infectious nbrs
          Rprintf("nInfInds %d\n", nInfInds);        
        }
      //
      //   // Infect nbr
       infProb = beta * (1.0 * nInfInds) / (1.0 * (nNbrs + 1));
       // Rprintf("infProb %.2f\n", infProb);
        isNewInf = as<int>(rbinom(1, 1, infProb));
      if(isNewInf == 1){
           SMax[nbrSusInd] = tt-1;
        //Rf_PrintValue(SMax);
        susVec[nbrSusInd] = -1; // no longer susceptible
        infVec[nbrSusInd] = nbrSusInd; // tell model agent is infectious for next step
         }

      }
    }
    
  }
  U(0, _ ) = newA0;
  U(1, _ ) = SMax;
  U(2, _ ) = IMax;
  
  
  return U;
  
  
}





// LoopType is 0 for looping over infectious and 1 for over susceptible (but smartly, looking at who was infected)
// Return the list of U arrays for each iteration L
// [[Rcpp::export]]
List AMSIR(int L, int T,
           IntegerVector A0,
           List nbrList,
           double beta,
           double gamma){
  
  List UList(L);
  for(int ll=0; ll < L; ll++){
    if((ll + 1) % 10 == 0){
      Rprintf("Iteration %d\n", ll + 1);
    }
  
      UList[ll] = AMSIR_sus_inner(ll, T,
                                  A0, nbrList,
                                  beta, gamma);
      
  }
  return UList;
  
} 


// 
// 
// /*** R
// 
// N <- 1000
// a0 <- rep(0, N)
// a0[25:75] <- 1
// T <- 100
// beta <- .1
// gamma <- .03
// nbr_list <- lapply(1:N, function(ii) (1:N)[-ii])
// L <- 100
// 
// t <- proc.time()[3]
// out <-  AMSIR(L, T, a0,
//       nbr_list, beta, gamma)
// proc.time()[3] - t
// 
// */
