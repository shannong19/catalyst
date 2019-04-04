## SKG
## March 26, 2019
## Plot object from am-sir.R


#' Plot the output of am_sir()
#'
#' @param X X output of result of am_sir()
#' @param cols vector of colors
#' @param title title to give ggplot
#' @param subtitle subtitle to give ggplot
#' @param obs data frame of observations with column names (t, S, I, R) .  Default is NULL, and so no observations are plotted
#' @param plot_var logical.  Do we plot CI?  Default is TRUE
#' @param labs model type labels.  Default is NULL
#' @param color_guide_name  Default is "Model Avg."
#' @return ggplot and the plotting df
am_plot_mean_var <- function(X, obs = NULL,
                             cols = c("blue", "red", "yellow3"),
                             title = "Hagelloch AM Simulations",
                             subtitle = "",
                             plot_var = TRUE,
                             labs = NULL,
                             color_guide_name = "Model Avg."){
    if(!("model"  %in% colnames(X))){
        X <- as.data.frame(X)
        X$model <- 1
        n_models <- 1
    } else{
        n_models <- length(unique(X$model))
    }

    if(is.null(labs)){
        labs <- 1:length(cols)
    }
    N <- X$S[1] + X$I[1] + X$R[1]
    Xm <- reshape2::melt(X, id.vars = c("t", "ll", "model"))
    ## Extract mean and variance plyr 4eva
    df <- plyr::ddply(Xm, .var = c("t", "variable", "model"),
                      .fun = function(df){
                          c("mean" = mean(df$value),
                            "var" = var(df$value))
                      })
    g <- ggplot2::ggplot() +
        ggplot2::geom_line(data = df,
                         ggplot2::aes(x = t, y = mean / N * 100, group = factor(model),
                                      col = factor(model)),
                         size = 1) +
        ggplot2::scale_color_manual(values = cols[1:n_models],
                                    name = color_guide_name, labels = labs) +
        ggplot2::scale_fill_manual(values = cols[1:n_models],  guide = FALSE) +
        my_theme() +
        ggplot2::facet_wrap(~variable, nrow = 3)
    if(plot_var){
        g <- g + 
            ggplot2::geom_ribbon(data = df, ggplot2::aes(x = t, group = factor(model),
                                                         ymin = mean / N * 100 - 2 * sqrt(var)/ N * 100,
                                                         ymax = mean / N * 100 + 2 * sqrt(var)/ N * 100,
                                                         fill = factor(model)),
                                 col = NA, alpha = .3)
    }
    g <- g + 
        ggplot2::labs(x = "Time",
                      y = "% of individuals",
                      title = title,
                      subtitle = subtitle,
                      color = "Model Avg.")

    if(!is.null(obs)){
        cn <- c("t", "S", "I", "R")
        obs_m <- reshape2::melt(obs[, cn], id.vars = "t")
        g <- g + ggplot2::geom_point(data = obs_m,
                                     ggplot2::aes(x = t, y = value / N * 100, group = variable))

    }

    
    print(g)
    return(list(g = g, plot_df = df))

}
