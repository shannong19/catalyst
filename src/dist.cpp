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
  List nbrList(n_obs);
  for(int ii=0; ii < n_obs; ii++){
    if(ii % 10000 == 0){
      Rprintf("agent %d\n", ii);
    }
    counter = 0;
    for(int jj=0; jj < n_obs; jj++){
      if(ii != jj){
        xlon = M(ii, 0);
        xlat = M(ii,1);
        ylon = M(jj, 0);
        ylat = M(jj, 1);
        dist = dist_haversine_rcpp(xlon, xlat, ylon, ylat);
      //  Rprintf("ii %d; jj %d; dist %.2f\n", ii, jj, dist);
        if(dist <= thresh){
          nbrVec[counter] = jj;
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





