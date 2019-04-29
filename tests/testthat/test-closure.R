context("testing closure routines")

test_that("removeClosureNbrs", {

    nbr_list <- list(c(1:4),
                     c(0, 2, 3),
                     c(0:1, 3:4),
                     c(4),
                     c(0:3))
    cat_inds <- 0:4

    out <- removeClosureNbrs(nbr_list, cat_inds)
    expect_equal(length(out), length(nbr_list))
    out_exp <- list(-1, -1, -1, -1, -1)
    expect_equal(out_exp, out)


    ## changing cat inds
    nbr_list <- list(c(1:4),
                     c(0, 2, 3),
                     c(0:1, 3:4),
                     c(4),
                     c(0:3))
    cat_inds <- 0
    out <- removeClosureNbrs(nbr_list, cat_inds)
    out_exp <- c(1:4)
    expect_equal(out_exp, out[[1]])

})


test_that("addClosureNbrs", {

    nbr_list <- list(-1, -1, -1, -1, -1)
    cat_inds <- 0:4

    out <- addClosureNbrs(nbr_list, cat_inds)
    expect_equal(length(out), length(nbr_list))
    out_exp <- list(c(1:4),
                    c(0, 2:4),
                    c(0:1, 3:4),
                    c(0:2, 4),
                    c(0:3))
    expect_equal(out_exp, out)


    ## changing cat inds
    nbr_list <- list(c(1:4),
                     c(0, 2, 3),
                     c(0:1, 3:4),
                     c(4),
                     c(0:3))
    cat_inds <- 0
    out <- addClosureNbrs(nbr_list, cat_inds)
    expect_equal(out, nbr_list)

    ## changing cat inds
    nbr_list2 <- list(c(2:4),
                     c(2, 3),
                     c(0:1, 3:4),
                     c(4),
                     c(0:3))
    cat_inds <- c(0,1)
    out <- addClosureNbrs(nbr_list2, cat_inds)
    expect_equal(out, nbr_list)

    
    ## changing cat inds
    nbr_list2 <- list(c(-1),
                     c(-1),
                     c(0:1, 3:4),
                     c(4),
                     c(0:3))
    cat_inds <- c(0,1)
    out_exp <- list(c(1),
                    c(0),
                    nbr_list2[[3]],
                    nbr_list2[[4]],
                    nbr_list2[[5]])
    out <- addClosureNbrs(nbr_list2, cat_inds)
    expect_equal(out, out_exp)

})


test_that("update_closure", {

    ## Set up
    nbr_list <- list(c(1:4),
                     c(0, 2, 3),
                     c(0:1, 3:4),
                     c(4),
                     c(0:3))
    inf_inds <- c(0)
    closure_thresh <- list(c(0, .2))
    closure_vars <- "school"
    closure_inds_list <- list(
        list(c(0, 1),
             c(2, 3, 4)
             ))
    closure_time <- list(c(0, 0))
    closure_max_T <- list(c(5, 5))
    out <- update_closure(nbr_list,
                          inf_inds,
                          closure_thresh,
                          closure_vars,
                          closure_inds_list,
                          closure_time,
                          closure_max_T)
    expect_equal(out$new_closure_time[[1]], c(1,0))

    ## Set up
    nbr_list <- list(c(1:4),
                     c(0, 2, 3),
                     c(0:1, 3:4),
                     c(4),
                     c(0:3))
    inf_inds <- c(0, 2)
    closure_thresh <- list(c(0, 1/3))
    closure_vars <- "school"
    closure_inds_list <- list(
        list(c(0, 1),
             c(2, 3, 4)
             ))
    closure_time <- list(c(0, 0))
    closure_max_T <- list(c(5, 5))
    out <- update_closure(nbr_list,
                          inf_inds,
                          closure_thresh,
                          closure_vars,
                          closure_inds_list,
                          closure_time,
                          closure_max_T)
    expect_equal(out$new_closure_time[[1]], c(1,0))


       ## Set up
    nbr_list <- list(c(1:4),
                     c(0, 2, 3),
                     c(0:1, 3:4),
                     c(4),
                     c(0:3))
    inf_inds <- c(0, 2)
    closure_thresh <- list(c(0, .1))
    closure_vars <- "school"
    closure_inds_list <- list(
        list(c(0, 1),
             c(2, 3, 4)
             ))
    closure_time <- list(c(0, 0))
    closure_max_T <- list(c(5, 5))
    out <- update_closure(nbr_list,
                          inf_inds,
                          closure_thresh,
                          closure_vars,
                          closure_inds_list,
                          closure_time,
                          closure_max_T)
    expect_equal(out$new_closure_time[[1]], c(1,1))


    ## Set up
    nbr_list <- list(c(1:4),
                     c(0, 2, 3),
                     c(0:1, 3:4),
                     c(4),
                     c(0:3))
    inf_inds <- c(0, 2)
    closure_thresh <- list(c(0, .1))
    closure_vars <- "school"
    closure_inds_list <- list(
        list(c(0, 1),
             c(2, 3, 4)
             ))
    closure_time <- list(c(1, 1))
    closure_max_T <- list(c(0, 0))
    out <- update_closure(nbr_list,
                          inf_inds,
                          closure_thresh,
                          closure_vars,
                          closure_inds_list,
                          closure_time,
                          closure_max_T)
    expect_equal(out$new_closure_time[[1]], c(1,1))
    exp_out <- nbr_list
    exp_out[[4]] <- c(2,4)
    expect_equal(out$new_nbrs, exp_out)

})


test_that("am_sir", {
## ugh

    })
