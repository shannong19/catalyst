
#' Calculate a statistic to summarize distance between two sets of points
#'
#' @param JT matrix of coordinates (the larger set)
#' @param J0 matrix of coordinates (the reference set) subset of JT
#' @param dist a function taking xlon, xlat, ylon, ylat
#' @param vector of minimum distances
mstat <- function(JT, J0, dist = dist_haversine_rcpp){
    JT <- as.matrix(JT)
    J0 <- as.matrix(J0)
    min_dists <- apply(JT, 1, mstat_inner, dist = dist, J0 = J0)

    return(min_dists)

}

#' Inner workings to return min dist
#' 
#' @param x row of mat with lon lat
#' @param J0 matrix of coordinates lon lat
#' @param dist a function
#' @return minimum distance between x and J0
mstat_inner <- function(x, J0, dist = dist){
  #print(x)
    dists <- apply(J0, 1, function(y){
        xlon = x[1]
        xlat = x[2]
        ylon = y[1]
        ylat = y[2]
        delta <- dist(xlon, xlat, ylon, ylat)
        return(delta)
    })
    return(min(dists))
}
