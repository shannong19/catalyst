## SKG
## Plotting SIR functions and themes
## March 6, 2019 Ash Wednesday.  Truly this is penance


#' Plot the t, X simulations
#'
#' @param catalyst_sims output from am_sir()
#' @param K number of total states
#' @param col_pal color palette of length K
#' @param states state names
#' @param xlab label for x
#' @param ylab label for y
#' @param title title
#' @param beta between 0 and 1
#' @param between 0 and  1
#' @param N total number of agents
#' @param L total number of runs
#' @param T total time steps
#' @param pretty logical.  Default is TRUE
#' @return ggplot graph
plot_sims <- function(catalyst_sims, K = 3,
                      col_pal = c("blue", "red", "yellow3"),
                      states = c("S", "I", "R"),
                      xlab = "Time",
                      ylab = "% of Individuals",
                      title = "SIR Simulations",
                      beta = .1,
                      gamma = .03,
                      N = 1000,
                      L = 100,
                      T = 100,
                      pretty = TRUE
                      ){

    df <- as.data.frame(sims$X_mat)
    dfm <- reshape2::melt(df, id.vars = c("t", "ll"))
    
    g <- ggplot2::ggplot(dfm, ggplot2::aes(x = t, y = value,
                                           group = paste(ll, variable),
                                           col = variable)) +
        ggplot2::geom_line() +
        ggplot2::scale_color_manual(values = col_pal, labels = states)
    g
                                                  
    
    

    return(g)
}
                      
