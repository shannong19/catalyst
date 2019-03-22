## SKG
## March 17, 2019 Happy St. Patrick's Day!!
## Loglinear model and plot.  Just using obs



#' Model log linear SIR
#'
#' @param obs data frame with t, S, I, and R,
#' @param CI_level between 0 and 1.  Default is .95
#' @param title title for ggplot
#' @param subtitle subtitle to be passed to ggplot
#' @param weighted_obs logical.  Default is TRUE
#' @param do_plot whether to plot the model or not
#' @return list with final model, the variance weights (if applicable), and plot (if applicable)
r0_loglinear_sir <- function(obs,
                             CI_level = .95,
                             title = "Log-linear regression line",
                             subtitle = "Measles in Hagelloch, Germany 1861",
                             weighted_obs = TRUE,
                             do_plot = TRUE){

    min_t <- which.min(df$t)[1]
    S0 <- df$S[min_t]
    N <- df$S[min_t] + df$I[min_t] + df$R[min_t]
    df$x <- df$R / N
    df$y <- -log(df$S / S0)
    ## Model 0, unweighted regression
    mod0 <- lm(y~x-1, data = df)
    ##
    out_mod <- mod0
    if(weighted_obs){
        res_df <- data.frame(res = mod0$res, x = df$x)
        weighted_mod <- lm(res^2~x-1, data = res_df)
        df$var <- predict(weight_mod)
        ## Weighted model (excluding Rt = 0 points)
        new_df <- df[df$x > 0,]
        mod1 <-lm(y~x-1, weights = 1 / var, data = new_df)
        out_mod <- mod1
    }
    pi <- predict(out_mod, interval = "prediction",
                  level = CI_level, newdata = data.frame(x = new_df$x),
                  weights = 1 / new_df$var)

    if(do_plot){
        est_df <- as.data.frame(pi)
        est_df$x <- new_df$x

        g <- ggplot2::ggplot() +
            ggplot2::geom_line(data = est_df,
                               ggplot2::aes(x = x, y = fit), col = "blue",
                               size = 2) +
            ggplot2::geom_ribbon(data = est_df,
                                 ggplot2::aes(x = x, ymin = lwr, max = upr),
                                 col = NA, fill = "blue", alpha = .3) +
            ggplot2::geom_point(data = df,
                                ggplot2::aes(x = x, y = y), size = 2)


        g <- g + my_theme() +
            ggplot2::labs(x = latex2exp::TeX("$R_t/N$"),
                          y = latex2exp::TeX("$-\\log (S_t / S_0 )$"),
                          title = title,
                          subtitle = subtitle)
        print(g)
                               
        

    } else{
        g <- NA
    }
    return(list(final_model = out_mod,
                CI_level = CI_level,
                g = g))

    
    


}
