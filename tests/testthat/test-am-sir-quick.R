context("quick")


test_that("findIfSus", {
    T <- 4
    A0 <- c(0, 1, 0, 0, 0, 1)
    IMax <- c(0, 1, 3, 3, 2, 3)
    ii <- 1
    out <- findIfSus(A0[ii], IMax[ii], T)
    expect_equal(out, 0)
    ##
    ii <- 3
    out <- findIfSus(A0[ii], IMax[ii], T)
    expect_equal(out, 1)



})


test_that("whichState", {
    A0 <- c(0, 1, 0, 0, 0, 1)
    out <- whichState(A0, 1)
    expect_equal(out, c(1, 5))
    ##
    A0 <- c(0, 0, 0, 0, 0, 0)
    out <- whichState(A0, 1)
    expect_equal(out, -1)


})

test_that("AMSIR",{

    

    a0 <- c(0, 1, 1)
    L <- 1
    T <- 5
    nbr_list = list(c(1, 2),
                    c(0, 2),
                    c(0, 1))
    beta <- 0
    gamma <- 1
    loopType <- 0
    out <- AMSIR(L, T, a0,
                 nbr_list, beta, gamma,
                 loopType)
    expect_equal(class(out), "list")
    expect_equal(length(out), 1)    
    exp_U <- matrix(c(a0,
                      T-1, T-1, T-1,
                      T-1, 0, 0), byrow = TRUE, nrow = 3)
    expect_equal(exp_U, out[[1]])
    ##
    a0 <- c(0, 1, 1)
    L <- 10
    T <- 5
    nbr_list = list(c(1, 2),
                    c(0, 2),
                    c(0, 1))
    beta <- 0
    gamma <- 1
    loopType <- 0
    out <- AMSIR(L, T, a0,
                 nbr_list, beta, gamma,
                 loopType)
    expect_equal(class(out), "list")
    expect_equal(length(out), L)
    expect_equal(exp_U, out[[L]])

 

}
)



test_that("quick am sir",{

    #devtools::load_all("~/catalyst")
    

    a0 <- c(0, 1, 1)
    ll <- 1
    T <- 5
    nbr_list = list(c(1, 2),
                    c(0, 2),
                    c(0, 1))
    beta <- 0
    gamma <- 1

    out <- AMSIR_inf_inner(ll, T, a0,
                           nbr_list, beta, gamma)
    exp_U <- matrix(c(a0,
                      T-1, T-1, T-1,
                      T-1, 0, 0), byrow = TRUE, nrow = 3)
    expect_equal(exp_U, out)

    ## test the quick and dirty AM, looping over infectious

    ## where inf has no neighbors
    a0 <- c(1, 0, 0)
    ll <- 1
    T <- 5
    nbr_list = list(c(-1),
                    c(2),
                    c(1))
    beta <- 1
    gamma <- 1
    out <- AMSIR_inf_inner(ll, T, a0,
                           nbr_list, beta, gamma)
    exp_U <- matrix(c(a0,
                      T-1, T-1, T-1,
                      0, T-1, T-1), byrow = TRUE, nrow = 3)
    expect_equal(exp_U, out)

}
)



test_that("quick am sir loop over sus",{

    #devtools::load_all("~/catalyst")
    

    a0 <- c(0, 1, 1)
    ll <- 1
    T <- 5
    nbr_list = list(c(1, 2),
                    c(0, 2),
                    c(0, 1))
    beta <- 0
    gamma <- 1

    out <- AMSIR_sus_inner(ll, T, a0,
                           nbr_list, beta, gamma)
    exp_U <- matrix(c(a0,
                      T-1, T-1, T-1,
                      T-1, 0, 0), byrow = TRUE, nrow = 3)
    expect_equal(exp_U, out)

    ## test the quick and dirty AM, looping over infectious

    ## where inf has no neighbors
    a0 <- c(1, 0, 0)
    ll <- 1
    T <- 5
    nbr_list = list(c(-1),
                    c(2),
                    c(1))
    beta <- 1
    gamma <- 1
    out <- AMSIR_sus_inner(ll, T, a0,
                           nbr_list, beta, gamma)
    exp_U <- matrix(c(a0,
                      T-1, T-1, T-1,
                      0, T-1, T-1), byrow = TRUE, nrow = 3)
    expect_equal(exp_U, out)

}
)



## Ye old favorite SIR
test_that("favorite SIR", {


    devtools::load_all("~/catalyst")
    ## Really slow cuz there are so many neighbors
    N <- 1000
    a0 <- rep(0, N)
    a0[1:50] <- 1
    T <- 100
    beta <- .1
    gamma <- .03
    nbr_list <- lapply(1:N, function(ii) (1:N)[-ii])
    L <- 100
    out <-  AMSIR(L, T, a0,
                 nbr_list, beta, gamma)

    test <- do.call('rbind',
                    lapply(1:length(out), function(ii){
                        mat <- out[[ii]]
                        df <- as.data.frame(UtoX_SIR(mat, T))
                        colnames(df) <- c("S", "I", "R")
                        df$t <- 0:(T-1)
                        df$ll <- ii
                        df <- df[, c("t", "S", "I", "R", "ll")]
                        return(df)
                    })
                    )
    df <- plyr::ddply(.data = test, .var = c("t"),
                      .fun = function(df){
                          S_mean = mean(df$S)
                          I_mean = mean(df$I)
                          R_mean = mean(df$R)
                          out <- data.frame(S_mean = S_mean,
                                            I_mean = I_mean,
                                            R_mean = R_mean)
                          return(out)
                      })

    plot(df$t, df$S_mean, type = "l", col = "blue", ylim = c(0, 1000))
    lines(df$t, df$I_mean, col = "red")
    lines(df$t, df$R_mean, col = "yellow3")
    
})


test_that("countIntersect", {

    devtools::load_all("~/catalyst")

    infVec <- c(-1, 1, -1, 3)
    nbrOfSusInsd <- c(1, 2)
    
    

})
