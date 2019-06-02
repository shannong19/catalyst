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


mstat <- function(JT, J0, dist = dist_haversine_rcpp){
    JT <- as.matrix(JT)
    J0 <- as.matrix(J0)
    min_dists <- apply(JT, 1, mstat_inner, J0 = J0,
                       dist = dist)

}


mstat_inner <- function(x, J0, dist = dist_haversine_rcpp){
    dists <- apply(J0, 1, function(y){
        xlon <- x[1]
        xlat <- x[2]
        ylon <- y[1]
        ylat <- y[2]
        dist <- dist(xlon, xlat, ylon, ylat)
        return(dist)

    })
    return(min(dists))
}
