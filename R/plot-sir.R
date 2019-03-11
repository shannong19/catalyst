## SKG
## Plotting SIR functions and themes
## March 6, 2019 Ash Wednesday.  Truly this is penance


#' Plot the observed
#'
#' @param obs data frame with t, S, I, and R columns (as percents out of total population N)
#' @param ests data frame with t, S_mean, I_mean, and R_mean columns along with est_type, and optionally S_var, I_var, R_var (as percents out of total population N)
#' @param plot_type "state" for regular SIR vs. t (faceted), "loglinear" for loglinear, and "ternary" for ternary
#' @param CI logical.  Default is FALSE.  Should we plot confidence intervals?
#' @param model_cols color palette
#' @return a ggplot
plot_ests <- function(obs, ests, plot_type = "state",
                      CI = FALSE,
                      model_cols = c("black", "blue", "red"),
                      xlab = "Time",
                      ylab = "% of Individuals",
                      title = "Observed data and fitted models",
                      CI_lab = NULL,
                      data_type = "Type",
                      pretty = TRUE,
                      model_names = c("Observed", "Model 1", "Model 2")
                      ){

    ## format obs
    obs_df <- format_obs(obs, plot_type, CI)
    ests_df <- format_ests(ests, plot_type, CI)
    df <- rbind(obs_df, ests_df)

    if(plot_type == "state"){
        g <- plot_ests.state(df)
    } else if(plot_type == "loglinear"){
        g <- plot_ests.loglinear(df)
    } else if(plot_type == "ternary"){
        g <- plot_ests.ternary(df)
    } else{
        stop("Select one of 'state', 'loglinear', or ternary'")
    }
    if(pretty){
        g <-  g + my_theme()
    }
    if(!is.null(CI_lab)){
        title <-  paste(title, "with", CI_lab)
    }

    g <- g + ggplot2::labs(x = xlab, y = ylab,
                  title = title) +
        ggplot2::scale_colour_manual(name = data_type,
                                       labels = model_names,
                                       values = model_cols) +
        ggplot2::scale_linetype_discrete(name = data_type,
                                       labels = model_names,) +
        ggplot2::scale_fill_manual(name = data_type,
                                     labels = model_names,
                                     values = model_cols)
    print(g)
        
    return(g)
}

#' Format observation data frame for plotting function
#'
#' @param ests data frame with t, S_mean, I_mean, and R_mean columns along with est_type, and optionally S_var, I_var, R_var (as percents out of total population N) and data_type
#' @param plot_type "state" for regular SIR vs. t (faceted), "loglinear" for loglinear, and "ternary" for ternary
#' @param CI logical.  Default is FALSE.  Should we plot confidence intervals?
#' @return data frame with columns "t", "obs", "state" (one of "S", "I", or "R"), "mean", and optionally "var",  along with "data_type"
format_ests <- function(ests, plot_type, CI){

    if(plot_type == "state"){
        var_order <- c("t", "obs", "state", "mean")
        if(CI){
            var_order <- c(var_order, "var", "data_type")
        } else{
            var_order <- c(var_order, "data_type")
        }


        df <- ests

        df$obs <- NA

        df$data_type <- factor(df$data_type,
                               labels = 1:length(unique(df$data_type)))

        ## Extract the mean and melt
        mean_vars <- grep("mean", colnames(df), value = TRUE)
        mean_df <- df[, c("t", "data_type", "obs", mean_vars)]
        mean_df_melt <- reshape2::melt(mean_df,
                                       id.vars = c("t", "data_type",
                                                   "obs"),
                                       value.name = "mean")
        mean_df_melt$state <- gsub("_mean", "", mean_df_melt$variable)
        ## Extract the variance and melt
        if(CI){
            var_vars <- grep("var", colnames(df), value = TRUE)
            var_df <- df[, c("t", "data_type", "obs", var_vars)]
            var_df_melt <- reshape2::melt(var_df,
                                           id.vars = c("t", "data_type",
                                                       "obs"),
                                           value.name = "var")
            var_df_melt$state <- gsub("_var", "", var_df_melt$variable)
            df <- mean_df_melt
            df$var <- var_df_melt$var
            df <- df[, var_order]
        } else{
            df <- mean_df_melt[, var_order]
        }
        

    } else if(plot_type == "loglinear"){

    } else if(plot_type == "ternary"){

    }


    df <- df[, var_order]
    return(df)
    

}



#' Format observation data frame for plotting function
#'
#' @param obs data frame with t, S, I, and R columns
#' @param plot_type "state" for regular SIR vs. t (faceted), "loglinear" for loglinear, and "ternary" for ternary
#' @param CI logical.  Default is FALSE.  Should we plot confidence intervals?
#' @return data frame with columns "t", "obs", "state" (one of "S", "I", or "R"), "mean", and optionally "var",  along with "data_type"
format_obs <- function(obs, plot_type, CI){
    if(plot_type == "state"){
        var_order <- c("t", "obs", "state", "mean")
        if(CI){
            var_order <- c(var_order, "var", "data_type")
        } else{
            var_order <- c(var_order, "data_type")
        }

        df <- reshape2::melt(obs, id.vars = c("t"))
        colnames(df) <- c("t", "state", "obs")
        df$mean <- NA
        if(CI){
            df$var <- NA
        }
        df$data_type <- 0
    } else if(plot_type == "loglinear"){

    } else if(plot_type == "ternary"){

    }

    df <- df[, var_order]


    return(df)
    

}

#' Plot the ggplot estimates for the regular states
#' 
#' @param df data frame with columns "t", "obs", "state" (one of "S", "I", or "R"), "mean", and optionally "var",  along with "data_type"
#' @param pretty logical.  Default is TRUE
#' @return ggplot
plot_ests.state <- function(df, pretty = TRUE){
    g <- ggplot2::ggplot(data = df,
                         ggplot2::aes(x= t, group = data_type)) +
        ggplot2::facet_wrap(~state, nrow = 3) +
        ggplot2:: geom_point(ggplot2::aes(y = obs,
                                          col = data_type),
                             col = "black", size = 1)
    if("var" %in% colnames(df)){
        g <- g +
            ggplot2::geom_ribbon(
                         ggplot2::aes(ymin = mean - 2 * sqrt(abs(var)),
                                      ymax = mean + 2 * sqrt(abs(var)),
                                      fill = data_type),
                         col = NA, alpha = .2)
        
    }
    g <- g + ggplot2::geom_line(ggplot2::aes(y = mean,
                                             col = data_type,
                                             linetype = data_type))

    return(g)
    
        
}

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
#' @export
plot_sims.reg <- function(catalyst_sims, K = 3,
                      col_pal = c("blue", "red", "yellow3"),
                      states = c("S", "I", "R"),
                      xlab = "Time",
                      ylab = "% of Individuals",
                      title = "SIR Simulations",
                      beta = .1,
                      gamma = .03,
                      N = 1000,
                      S0 = 950,
                      I0 = 50,
                      L = 100,
                      T = 100,
                      pretty = TRUE
                      ){

    df <- as.data.frame(sims$X_mat)
    dfm <- reshape2::melt(df, id.vars = c("t", "ll"))
    
    g <- ggplot2::ggplot(dfm, ggplot2::aes(x = t, y = value / N * 100,
                                           group = paste(ll, variable),
                                           col = variable)) +
        ggplot2::geom_line() +
        ggplot2::scale_color_manual(values = col_pal, labels = states,
                                    name = "State") +
        ggplot2::labs(title = title,
                      subtitle = latex2exp::TeX(sprintf("N = %d, $S_0$ = %d, $I_0$ = %d, T = %d, L = %d, $\\beta = %.2f$, $\\gamma = %.2f$", N, S0, I0, T, L, beta, gamma)),
                      x = xlab,
                      y = ylab
                      )

    if(pretty)
        g <- g + my_theme()
    
    print(g)
                                                  
    
    return(g)
}
                      

#' My custom ggplot background
#'
#' @return ggplot theme
my_theme <- function(){
    ggplot2::theme_bw() +
        ggplot2::theme(
                     axis.text.x = ggplot2::element_text(size = 16,
                                                         family = "Palatino"),
                     axis.text.y= ggplot2::element_text(size = 16,
                                                        family = "Palatino"),
                     axis.title.x= ggplot2::element_text(size = 18,
                                                         family = "Palatino"),
                     axis.title.y= ggplot2::element_text(size = 18,
                                                         family = "Palatino"),
                     plot.title = ggplot2::element_text(size = 24,
                                                        family = "Palatino"),
                     legend.title = ggplot2::element_text(size = 20,
                                                          family = "Palatino"),
                     legend.text = ggplot2::element_text(family = "Palatino",
                                                         size = 16),
                     legend.key.size = ggplot2::unit(3, "line"),
                     plot.subtitle = ggplot2::element_text(size=16,
                                                           family = "Palatino"),
                     strip.text = ggplot2::element_text(size = 16,
                                                        family = "Palatino")
                 ) +
            ggplot2::theme(legend.position = "bottom") 
      
}


#' Turn raw estimates to percent of N
#'
#' @param df with columns t, S, I, and R
#' @param N total number of individuals
#' @return transformed df
df_to_pct <- function(df, N){
    mean_vars <- grep("mean", colnames(df), value = TRUE)
    mean_vars <- unique(c(mean_vars,
                          colnames(df)[which(colnames(df) %in% c("S", "I", "R"))]
                          ))
    var_vars <- grep("var", colnames(df), value = TRUE)
    df[, mean_vars] <- df[, mean_vars] * 100 / N
    df[, var_vars] <- df[, var_vars] * (100 / N)^2

    return(df)

}

#' Plot the mean and variance of simulations
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
plot_mean_var.reg <- function(catalyst_sims, K = 3,
                          col_pal = c("blue", "red", "yellow3"),
                          states = c("S", "I", "R"),
                          xlab = "Time",
                          ylab = "% of Individuals",
                          title = "SIR Simulations",
                          beta = .1,
                          gamma = .03,
                          N = 1000,
                          S0 = 950,
                          I0 = 50,
                          L = 100,
                          T = 100,
                          pretty = TRUE
                          ){

    df <- as.data.frame(catalyst_sims$x)

    dfm <- reshape2::melt(df, id.vars = c("t", "ll"))
    summary_df <- plyr::ddply(dfm, .var = c("t", "variable"),
                            .fun  = function(df){
                                c(mean = mean(df$value), var = var(df$value))
                            })


    g <- ggplot2::ggplot(summary_df, ggplot2::aes(x = t, y = mean / N * 100,
                                           group = variable,
                                           col = variable)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = mean / N * 100 - 2 * sqrt(var) / N * 100,
                                 ymax = mean / N * 100 + 2 * sqrt(var) / N * 100,
                                 fill = variable), alpha = .7, col = NA) +
        ggplot2::geom_line(size = 2) +
        ggplot2::scale_color_manual(values = col_pal, labels = states,
                                    name = "State") +
        ggplot2::scale_fill_manual(values = col_pal, name = "State") +
        ggplot2::labs(title = title,
                      subtitle = latex2exp::TeX(sprintf("N = %d, $S_0$ = %d, $I_0$ = %d, T = %d, L = %d, $\\beta = %.2f$, $\\gamma = %.2f$", N, S0, I0, T, L, beta, gamma)),
                      x = xlab,
                      y = ylab
                      )

    if(pretty)
        g <- g + my_theme()
    
    print(g)
    return(g)

}




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
#' @export
plot_sims.loglin <- function(catalyst_sims, K = 3,
                      col_pal = c("blue", "red", "yellow3"),
                      states = c("S", "I", "R"),
                      xlab = "Time",
                      ylab = "% of Individuals",
                      title = "SIR Simulations",
                      beta = .1,
                      gamma = .03,
                      N = 1000,
                      S0 = 950,
                      I0 = 50,
                      L = 100,
                      T = 100,
                      pretty = TRUE
                      ){

    df <- as.data.frame(sims$X_mat)
    df$x <- df$R / N
    df$y <- -log(df$S / df$S[1])
    dfm <- df# reshape2::melt(df, id.vars = c("t", "ll"))
    
    g <- ggplot2::ggplot(dfm, ggplot2::aes(x = x, y = y, group = ll, col = ll),
                         alpha = .4) +
        ggplot2::geom_line() +
        ggplot2::scale_color_gradient(guide = FALSE, low = "black", high = "white") + 
        ggplot2::labs(title = title,
                      subtitle = latex2exp::TeX(sprintf("N = %d, $S_0$ = %d, $I_0$ = %d, T = %d, L = %d, $\\beta = %.2f$, $\\gamma = %.2f$", N, S0, I0, T, L, beta, gamma)),
                      x = xlab,
                      y = ylab
                      )

    if(pretty)
        g <- g + my_theme()
    
    print(g)
                                                  
    
    return(g)
}
                      


#' Plot the mean and variance of simulations
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
plot_mean_var.reg <- function(catalyst_sims, K = 3,
                          col_pal = c("blue", "red", "yellow3"),
                          states = c("S", "I", "R"),
                          xlab = "Time",
                          ylab = "% of Individuals",
                          title = "SIR Simulations",
                          beta = .1,
                          gamma = .03,
                          N = 1000,
                          S0 = 950,
                          I0 = 50,
                          L = 100,
                          T = 100,
                          pretty = TRUE
                          ){

    df <- as.data.frame(sims$X_mat)
    df$x <- df$R / N
    df$y <- -log(df$S / df$S[1])
    dfm <- df
    summary_df <- plyr::ddply(dfm, .var = c("t", "variable"),
                              .fun  = function(df){
                                  c(mean = mean(df$value), var = var(df$value))
                              })


    g <- ggplot2::ggplot(summary_df, ggplot2::aes(x = t, y = mean / N * 100,
                                           group = variable,
                                           col = variable)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = mean / N * 100 - 2 * sqrt(var) / N * 100,
                                 ymax = mean / N * 100 + 2 * sqrt(var) / N * 100,
                                 fill = variable), alpha = .7, col = NA) +
        ggplot2::geom_line(size = 2) +
        ggplot2::scale_color_manual(values = col_pal, labels = states,
                                    name = "State") +
        ggplot2::scale_fill_manual(values = col_pal, name = "State") +
        ggplot2::labs(title = title,
                      subtitle = latex2exp::TeX(sprintf("N = %d, $S_0$ = %d, $I_0$ = %d, T = %d, L = %d, $\\beta = %.2f$, $\\gamma = %.2f$", N, S0, I0, T, L, beta, gamma)),
                      x = xlab,
                      y = ylab
                      )

    if(pretty)
        g <- g + my_theme()
    
    print(g)
    return(g)

}
