## SKG
## May 20, 2019
## Assign infections to covariate DF based on covariates


#' Assign infection and recovery date to agent df
#'
#' @param full_df full data frame with ALL agents of a population
#' @param inf_df data frame with infection and recovery date of all infected agents
#' @param assign_fxn which function we use to assign which agents were infected.  Default is by age
#' @param arg_list list of arguments to give to assign_fxn
#' @return 2 column data frame of I and R which are the infection date and recovery date, respectively.
assign_inf_recov_dates <- function(full_df,
                                   inf_df,
                                   assign_fxn = assign.byage,
                                   arg_list){

    new_df <- assign_fxn(full_df, inf_df, arg_list)
    return(new_df)

}


#' Assign infection and recovery date to agent df based on age
#'
#' @param full_df full data frame with ALL agents of a population
#' @param inf_df data frame with infection and recovery date of all infected agents
#' @param assign_fxn which function we use to assign which agents were infected.  Default is by age
#' @param arg_list list of arguments to give to assign_fxn.  entry is age_cutoff.  How close do we have to be to the reported age?
#' @return 2 column data frame of I and R which are the infection date and recovery date, respectively.
#' @details Requires full_df and inf_df to have an "AGE" column and inf_df must have an I and R column
assign.byage <- function(full_df,
                         inf_df,
                         arg_list = list(age_cutoff = 1,
                                         sequential = FALSE)){

    full_df$ind <- 1:nrow(full_df) # Need to maintain original order
    assigned_df <- plyr::ddply(.data = inf_df, .var = c("AGE"),
                               .fun = function(df){
                                   age <- unique(df$AGE)
                                   lower <- age - arg_list$age_cutoff
                                   upper <- age + arg_list$age_cutoff
                                   sub_inds <- which(full_df$Age >= lower |
                                                     full_df$Age <= upper)
                                   n <- nrow(df)
                                   sampled_inds <- sample(sub_inds, size = n)
                                   new_sub_df <- data.frame(ind = sampled_inds,
                                                            I = df$I,
                                                            R = df$R)
                                   return(new_sub_df)

                          })
    new_df <- dplyr::left_join(full_df, assigned_df, by = "ind")
    
    return(new_df)

}
