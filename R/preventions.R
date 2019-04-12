## SKG
## April 12, 2019
## CARNIVAL
## Saving children



#' Calculate list of neighbors from environment DATA frame
#'
#' @param env_df NxE data frame.  Default is NULL
#' @param vars variable names which we subset to when applying preventions.  eg isolation means only neighbors are household ID.
#' @return neighbors list if we were doing preventions for each agent
precompute_preventions_nbrs <- function(env_df, vars){
    if(!all(vars %in% colnames(env_df))) stop("These are not prevention variable options")
    df <- env_df[, vars]
    mat <- as.matrix(df)
    nbr_list <- initializeNeighbors(mat)
    return(nbr_list)


}


