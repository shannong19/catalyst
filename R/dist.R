## From http://www.guru-gis.net/fast-euclidean-distance-function/


#' Calculate Euclidean distance
#'
#' @param m Matrix where rows are observations and columns are features
#' @return matrix of distances
euc_dist <- function(m){
    mtm <- Matrix::tcrossprod(m)
    sq <- rowSums(m*m)
    sqrt(outer(sq,sq,"+") - 2*mtm)
} 
