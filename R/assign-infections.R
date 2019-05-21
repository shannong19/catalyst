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
                         arg_list = list(age_cutoff = 1)){

    stopifnot(nrow(full_df) >= nrow(inf_df)) ## Need to have agent set large enough
    full_df$ind <- 1:nrow(full_df) # Need to maintain original order
    inf_df$ebola_id <- 1:nrow(inf_df)

    ## Split by ages
    random_ages <- sample(unique(inf_df$AGE))
    sub_df <- full_df  ## Initially we sample indices from full data frame
    new_sub_list <- vector(mode = "list", length = length(unique(inf_df$AGE)))
    print(random_ages)
    for(ii in 1:length(random_ages)){
        ## Subset infection DF
        age <- random_ages[ii]
        if(is.na(age)){
            df <- inf_df[which(is.na(inf_df$AGE)),]
            sub_inds <- NULL
        } else {
            df <- inf_df[which(inf_df$AGE == age), ]
            lower <- age - arg_list$age_cutoff
            upper <- age + arg_list$age_cutoff
            sub_inds <- sub_df$ind[which(sub_df$AGE >= lower &
                                         sub_df$AGE <= upper)]
        }
        N <- nrow(df)
        
     #   print(paste("Age:", age, "; N", N))
       
        if(length(sub_inds) == 1){
            sampled_inds <- sub_inds
        }  else if(length(sub_inds) < N | is.na(age)){
            ## If there are not enough subjects to sample from, expand to whole (remaining) DF
            sub_inds <- sub_df$ind
            sampled_inds <- base::sample(x = sub_inds, size = N)
        } else {
            sampled_inds <- base::sample(x = sub_inds, size = N)
        }
        new_sub_list[[ii]] <- data.frame(ind = sampled_inds,
                                         I = df$I,
                                         R = df$R,
                                         AGE_inf = df$AGE,
                                         ebola_id = df$ebola_id)
        ## Get rid of sampled inds, cannot be picked again
        sub_df <- sub_df[ - which(sub_df$ind %in% sampled_inds), ]
    }
    assigned_df <- do.call('rbind', new_sub_list)
    print(length(unique(assigned_df$ind)))
    stopifnot(length(unique(assigned_df$ind)) == nrow(inf_df))
    new_df <- dplyr::left_join(full_df, assigned_df, by = "ind")
    new_df <- new_df[, - which(colnames(new_df) %in% "ind")]
    
    return(new_df)

}
