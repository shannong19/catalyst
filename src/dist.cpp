#include <Rcpp.h>
using namespace Rcpp;
// Thanks this guy: https://www.btskinner.io/rworkshop/modules/hp_rcpp.html

// preprocessor replacements
#define e_r 3961.0 // radius of earth in miles  
#define m_pi 3.141592653589793238462643383280 // pi




// convert degrees to radians
// [[Rcpp::export]]
double deg_to_rad_rcpp(double degree) {
  return(degree * m_pi / 180.0);
}

// compute Haversine distance between two points
// [[Rcpp::export]]
double dist_haversine_rcpp(double xlon,
                           double xlat,
                           double ylon,
                           double ylat) {
  
  // return 0 if same point
  if (xlon == ylon && xlat == xlon) return 0.0;
  
  // convert degrees to radians
  xlon = deg_to_rad_rcpp(xlon);
  xlat = deg_to_rad_rcpp(xlat);
  ylon = deg_to_rad_rcpp(ylon);
  ylat = deg_to_rad_rcpp(ylat);
  
  // haversine distance formula
  double d1 = sin((ylat - xlat) / 2.0);
  double d2 = sin((ylon - xlon) / 2.0);
  return 2.0 * e_r * asin(sqrt(d1*d1 + cos(xlat) * cos(ylat) * d2*d2));  
}



// [[Rcpp::export]]
Rcpp::IntegerVector sample_int(int n, int min, int max) {
  Rcpp::IntegerVector pool = Rcpp::seq(min, max);
  std::random_shuffle(pool.begin(), pool.end());
  return pool[Rcpp::Range(0, n - 1)];
}


// Get the neighbor list based on haversine distance
// 
// @param M matrix of size Nx2 where N is the number of agents
// @param thresh distance in miles of what is an acceptable neighbor
// @param max_nbrs maximum number of neighbors to consider (the space saver)
// [[Rcpp::export]]
List nbrsByDist(NumericMatrix M, double thresh, int max_nbrs){
  int n_obs = M.nrow();
  int n_col = M.ncol();
  
  double dist;
  int counter=0;
  double xlon;
  double xlat;
  double ylon;
  double ylat;
  IntegerVector nbrVec(max_nbrs);
  IntegerVector permutedInds(n_obs);
  int permutedInd;
  List nbrList(n_obs);
  for(int ii=0; ii < n_obs; ii++){
    if(ii % 10000 == 0){
      Rprintf("agent %d\n", ii);
    }
    counter = 0;
    // Need to permute the indices
    permutedInds = sample_int(n_obs, 0, n_obs-1);
    for(int jj=0; jj < n_obs; jj++){
      permutedInd = permutedInds[jj];
      if(ii != permutedInd){
        xlon = M(ii, 0);
        xlat = M(ii,1);
        ylon = M(permutedInd, 0);
        ylat = M(permutedInd, 1);
        dist = dist_haversine_rcpp(xlon, xlat, ylon, ylat);
      //  Rprintf("ii %d; jj %d; dist %.2f\n", ii, jj, dist);
        if(dist <= thresh){
          nbrVec[counter] = permutedInd;
          counter = counter + 1;
        }
        if(counter >= max_nbrs){  // Have too many neighbors
          break;
        }
      
      }
    }
   // Rprintf("counter %d\n", counter);
    if(counter == 0){
      nbrList[ii] = -1; // no neighbors
    } else {
      IntegerVector finalNbr(counter);
      for(int jj=0; jj < counter; jj++){ // shrink the neighbor vector
        finalNbr[jj] = nbrVec[jj];
      }
      nbrList[ii] = finalNbr;
    }
  
  }
  
  return nbrList;
}




/*** R
N <- 4
samp <- sample_int(4, 0, N -1)
 */
