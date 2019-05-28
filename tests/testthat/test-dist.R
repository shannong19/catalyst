context("testing Rcpp functions for haversine distance")


test_that("nbrsByDist", {


    devtools::load_all("~/catalyst")
    mat <- matrix(c(0, 1,
                    2, 2,
                    2, 1,
                    1, 1), byrow = TRUE, ncol = 2)


    thresh <- 70
    max_nbrs <- 10
    out <- nbrsByDist(mat, thresh, max_nbrs)
    hdist <- dist_haversine_rcpp(0, 1, 1, 1)
    expect_equal(length(out), 4)
    expect_equal(out[[1]], 3)
################################3

    
    mat <- matrix(c(0, 1,
                    2, 2,
                    2, 1,
                    1, 1,
                    5, 5), byrow = TRUE, ncol = 2)
    thresh <- 70
    max_nbrs <- 10
    out <- nbrsByDist(mat, thresh, max_nbrs)
    hdist <- dist_haversine_rcpp(0, 1, 1, 1)
    expect_equal(length(out), 5)
    expect_equal(out[[1]], 3)
    expect_equal(out[[5]], -1)


})
