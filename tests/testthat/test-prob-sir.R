context("testing probabilitiy functions")

test_that("Group probabilities", {
    T <- 4
    X <- matrix(c(4, 1, 1,
                  1, 5, 10,
                  4, 1, 2), nrow = 3)
    par <- c(1, 0, .1)
    N <- 100
    out <- KM_probG(par, X, inf_nbrs = NULL, N)
    expect_equal(dim(out), c(2, 3))
    



})
