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
    Xm <- Xm %>% dplyr::group_by(variable, t, model) %>%
        dplyr::summarize(mean = mean(value), var = var(value))
    g <- ggplot2::ggplot() +
        ggplot2::geom_line(data = Xm,
                         ggplot2::aes(x = t, y = mean / N * 100, group = factor(model),
                                      col = factor(model)),
                         size = 1) +
        ggplot2::scale_color_manual(values = cols[1:n_models],
                                    name = color_guide_name, labels = labs) +
        ggplot2::scale_fill_manual(values = cols[1:n_models],  guide = FALSE) +
        my_theme() +
        ggplot2::facet_wrap(~variable, nrow = 3, scales = "free_y")
    if(plot_var){
        g <- g + 
            ggplot2::geom_ribbon(data = Xm, ggplot2::aes(x = t, group = factor(model),
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


#' Plot epidemic summary of SIR AM
#'
#' @param sims_list list of the "X" output from am_sir() function
#' @param N number of total agents
#' @param title character passed on to ggplot
#' @param subtitle character passed on to ggplot
#' @param summary_fxn either "mean" or "median".  How we summarize the summary statistics
#' @param many_groups logical.  Default is FALSE.  If TRUE, we will do different colors on them.  Corresponds to the "Type" column
#' @return list with two ggplots and summarized, plottable df
plot_epidemic_summary <- function(sims_list, N,
                                  title,
                                  subtitle,
                                  summary_fxn = "mean",
                                  many_groups = FALSE,
                                  cols = "black",
                                  type_name = "Type",
                                  labs = NULL,
                                  xmin = 0,
                                  xmax = 80,
                                  y1min = -25,
                                  y1max = 75,
                                  y2min = -50,
                                  y2max = 150,
                                  y3min = -5,
                                  y3max = 120,
                                  lty = 0,
                                  show_labels = FALSE,
                                  alpha = .3,
                                  do_facet = FALSE,
                                  do_color_ramp = FALSE,
                                  pal = "Greys",
                                  legend_name = latex2exp::TeX("$\\rho$"),
                                  bg = FALSE){

    df <- summarize_epidemic(sims_list, N, summary_fxn)

    ## Plot max t vs max I
    plot_list <- plot_epidemic_df(df, N, many_groups = many_groups,
                                  cols = cols,
                                  title,
                                  subtitle,
                                  type_name,
                                  labs,
                                  xmin = xmin,
                                  xmax = xmax,
                                  y1min = y1min,
                                  y1max = y1max,
                                  y2min = y2min,
                                  y2max = y2max,
                                  y3min = y3min,
                                  y3max = y3max,
                                  lty = lty,
                                  show_labels = show_labels,
                                  alpha = alpha,
                                  do_facet = do_facet,
                                  do_color_ramp = do_color_ramp,
                                  pal = pal,
                                  legend_name = legend_name,
                                  bg = bg)
   
    return(list(g1 = plot_list$g1, g2 = plot_list$g2,
                g3 = plot_list$g3,
                df = df))

    
}

#' Summarize list of AM simulations
#' 
#' @param sims_list list of the "X" output from am_sir() function
#' @param N total number of agents
#' @return summary data frame with columns mean_* and var_* for ("max_I", "sum_I", and "max_t", "len_I") and 95% CIs for each of these vars around the mean
summarize_epidemic <- function(sims_list, N, summary_fxn = "mean"){

    if(summary_fxn == "median"){
        sum_fxn <- median
    } else{
        sum_fxn = mean
    }

    X <- do.call('rbind', sims_list)
    if("mult" %in% colnames(X)){
        X$model <- X$mult
        df <- X[, - which(colnames(X) == "mult")]
    } else {
        df <- X
    }
    if(!("color_ramp" %in% colnames(X))){
        df$color_ramp <- 1
    }
    if(!("Type" %in% colnames(X))){
        df$Type <- "Type"
    } else {
        if(class(df$Type) != "factor"){
            df$Type <- factor(df$Type)
        }
    }


    vars <- c("model", "ll", "color_ramp", "Type")
    vars2 <- c("model", "color_ramp", "Type")


    
    ## Make into a function
    sum_df <- plyr::ddply(df, .var = vars,
                          .fun = function(df){
                              c(max_I = max(df$I),
                                sum_I = max(df$R),
                                max_t = df$t[which.max(df$I)],
                                len_I = max(which(df$I >= 1))
                                    )
                          })
    ave_sum_df <- plyr::ddply(sum_df, .var = vars2,
                              .fun = function(df){
                                  c(mean_max_I = sum_fxn(df$max_I),
                                    var_max_I = var(df$max_I),
                                    mean_sum_I = sum_fxn(df$sum_I),
                                    var_sum_I = var(df$sum_I),
                                    mean_max_t = sum_fxn(df$max_t),
                                    var_max_t = var(df$max_t),
                                    mean_len_I = sum_fxn(df$len_I),
                                    var_len_I = var(df$len_I)
                                    )
                              }) 
    ave_sum_df$tmin <-  (ave_sum_df$mean_max_t -
                         2 * sqrt(ave_sum_df$var_max_t))
    ave_sum_df$tmax <-  (ave_sum_df$mean_max_t +
                         2 * sqrt(ave_sum_df$var_max_t))
    ave_sum_df$Imin <-  (ave_sum_df$mean_max_I / N * 100 -
                         2 * sqrt(ave_sum_df$var_max_I) / N * 100)
    ave_sum_df$Imax <-  (ave_sum_df$mean_max_I / N * 100 +
                         2 * sqrt(ave_sum_df$var_max_I) / N * 100)
    ave_sum_df$sumImin <-  (ave_sum_df$mean_sum_I / N * 100 -
                            2 * sqrt(ave_sum_df$var_sum_I) / N * 100)
    ave_sum_df$sumImax <-  (ave_sum_df$mean_sum_I / N * 100 +
                            2 * sqrt(ave_sum_df$var_sum_I) / N * 100)
    ave_sum_df$lenImin <-  (ave_sum_df$mean_len_I / N * 100 -
                            2 * sqrt(ave_sum_df$var_len_I) / N * 100)
    ave_sum_df$lenImax <-  (ave_sum_df$mean_len_I  +
                            2 * sqrt(ave_sum_df$var_len_I) )
    return(ave_sum_df)


}


plot_epidemic_df <- function(df, N, many_groups, cols = "black",
                             title,
                             subtitle,
                             type_name = "Type",
                             labs = NULL,
                             xmin = 0,
                             xmax = 80,
                             y1min = -50,
                             y1max = 100,
                             y2min = -25,
                             y2max = 150,
                             y3min = -5,
                             y3max = 120,
                             lty = 1,
                             show_labels = FALSE,
                             alpha = .3,
                             do_facet = FALSE,
                             do_color_ramp = FALSE,
                             pal = "Greys",
                             legend_name = latex2exp::TeX("$\\rho$"),
                             bg = FALSE
                             ){




  
    col_guide <- NULL
    if(!many_groups){
           df$Type <- "Model"
           col_guide <- FALSE
           if(!("color_ramp" %in%colnames(df))){
               df$color_ramp <- df$Type
           }
           } else{
               if(!("color_ramp" %in% colnames(df))){
                   df$color_ramp <- df$Type
               }
    }

    if(is.null(labs)){
        labs <- unique(df$Type)
    }
    if(do_color_ramp){
        df$color_ramp <- as.numeric(as.character(df$model))
    } else {
        df$color_ramp <- factor(df$Type)
    }
 
     g1 <- ggplot2::ggplot(data = df,
                          ggplot2::aes(x = mean_max_I / N * 100,
                                       y = mean_max_t,
                                       group = color_ramp,
                                       col = color_ramp,
##                                       shape = Type,
                                       x0 = mean_max_I / N * 100,
                                       y0 = mean_max_t,
                                       a = 2 * sqrt(var_max_I) / N * 100,
                                       b = 2 * sqrt(var_max_t),
                                       angle = 0,
                                       fill = color_ramp)) +
        ggplot2::geom_point(size = 2) + my_theme() +
        ggplot2::labs(x = "Peak % infectious over all days",
             y = "Day of peak % infectious",
             title = title,
             subtitle = subtitle) +
        ggplot2::guides(color = col_guide, shape = col_guide, fill = col_guide) +
        ggforce::geom_ellipse(alpha = alpha,  size  =1, linetype = lty) +
        ggplot2::ylim(y1min, y1max) +
        ggplot2::xlim(xmin, xmax)
    if(show_labels){
        g1 <-  g1 + ggplot2::geom_text(ggplot2::aes(label = model), nudge_x = -3 / N * 100,
                                       nudge_y = 5, show_legend = FALSE) 
    }
    if(do_facet){
        g1 <- g1 + ggplot2::facet_wrap(~Type)
    }

    if(do_color_ramp){
        g1 <-  g1 +  ggplot2::scale_color_distiller(legend_name,
                                                    palette = pal, direction  =1 ) +
            ggplot2::scale_fill_distiller(legend_name, palette = pal, direction = 1 ) 
            
    } else{
        g1 <- g1 +  ggplot2::scale_color_manual(values = cols, name = type_name, labels = labs) +
            ggplot2::scale_fill_manual(values = cols, name = type_name, labels = labs) 
    }

        
        ## ggplot2::geom_errorbar(ggplot2::aes(ymin = tmin,
        ##                                     ymax = tmax)) +
        ## ggplot2::geom_errorbarh(ggplot2::aes(xmin = Imin,
        ##                                      xmax = Imax)) +
       


    ## Plot max total I vs max I
    g2 <- ggplot2::ggplot(data = df,
                          ggplot2::aes(x = mean_max_I / N * 100,
                                       y = mean_sum_I / N * 100,
                                       group = color_ramp,
                                       col = color_ramp,
                                    #   shape = Type,
                                       x0 = mean_max_I / N * 100,
                                       y0 = mean_sum_I / N * 100,
                                       a = 2 * sqrt(var_max_I) / N * 100,
                                       b = 2 * sqrt(var_sum_I) / N * 100,
                                       angle = 0,
                                       fill = color_ramp))
    if(bg){
        g2 <- g2 +         ggplot2::geom_point(size = 3, col = "black")
    }
    g2 <-  g2 + 
        ggplot2::geom_point(size = 2) + my_theme() +
        ggplot2::labs(x = "Peak % infectious over all days",
             y = "Final size (%)",
             title = title,
             subtitle = subtitle) +
        ggplot2::ylim(y2min, y2max) +
        ggplot2::xlim(xmin, xmax) + 
        ## ggplot2::geom_errorbar(ggplot2::aes(ymin = sumImin,
        ##                                     ymax = sumImax)) +
        ## ggplot2::geom_errorbarh(ggplot2::aes(xmin = Imin,
        ##                                      xmax = Imax)) +
        ggplot2::guides(color = col_guide, shape = col_guide, fill = col_guide)

    if(bg){
           g2 <-  g2 + ggforce::geom_ellipse(alpha = alpha, size  =2, linetype = lty, col = "black")
    }
    g2 <- g2 + 
        ggforce::geom_ellipse(alpha = alpha, size  =1, linetype = lty)
    if(show_labels){
        g2 <- g2 + ggplot2::geom_text(ggplot2::aes(label = model), nudge_x = -3/ N * 100,
                                      nudge_y = 5)
    }
    if(do_facet){
        g2 <- g2 + ggplot2::facet_wrap(~Type)
    }

    if(do_color_ramp){
        g2 <-  g2 +  ggplot2::scale_color_distiller(legend_name, palette = pal, direction  = 1) +
            ggplot2::scale_fill_distiller(legend_name, palette = pal, direction  = 1) 
    } else{
        g2 <- g2 +  ggplot2::scale_color_manual(values = cols, name = type_name, labels = labs) + 
            ggplot2::scale_fill_manual(values = cols, name = type_name, labels = labs) 
    }
################################################################3
    ## Plot length I vs max I
    g3 <- ggplot2::ggplot(data = df,
                          ggplot2::aes(x = mean_max_I / N * 100,
                                       y = mean_len_I ,
                                       group = color_ramp,
                                       col = color_ramp,
                                    #   shape = Type,
                                       x0 = mean_max_I / N * 100,
                                       y0 = mean_len_I ,
                                       a = 2 * sqrt(var_max_I) / N * 100,
                                       b = 2 * sqrt(var_len_I),
                                       angle = 0,
                                       fill = color_ramp)) +
        ggplot2::geom_point(size = 2) + my_theme() +
        ggplot2::labs(x = "Peak % infectious over all days",
             y = "Infection duration",
             title = title,
             subtitle = subtitle) +
        ggplot2::ylim(y3min, y3max) +
        ggplot2::xlim(xmin, xmax) + 
        ## ggplot2::geom_errorbar(ggplot2::aes(ymin = sumImin,
        ##                                     ymax = sumImax)) +
        ## ggplot2::geom_errorbarh(ggplot2::aes(xmin = Imin,
        ##                                      xmax = Imax)) +
        ggplot2::guides(color = col_guide, shape = col_guide, fill = col_guide) +
        ggforce::geom_ellipse(alpha = alpha, size  =1, linetype = lty)
    if(show_labels){
        g3 <- g3 + ggplot2::geom_text(ggplot2::aes(label = model), nudge_x = -3/ N * 100,
                                      nudge_y = 5)
    }
    if(do_facet){
        g3 <- g3 + ggplot2::facet_wrap(~Type)
    }

    if(do_color_ramp){
        g3 <-  g3 +  ggplot2::scale_color_distiller(legend_name, palette = pal, direction  = 1) +
            ggplot2::scale_fill_distiller(legend_name, palette = pal, direction  = 1) 
    } else{
        g3 <- g3 +  ggplot2::scale_color_manual(values = cols, name = type_name, labels = labs) + 
            ggplot2::scale_fill_manual(values = cols, name = type_name, labels = labs) 
    }


    gridExtra::grid.arrange(g1, g2, g3, nrow = 3)
    return(list(g1 = g1, g2 = g2, g3 = g3))
}

