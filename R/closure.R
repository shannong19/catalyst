
#' Update the neighbors list due to closure of a environmental locus
#'
#' @param nbr_list list/map where the key is the is the index of the agent and values are its neigbhors or -1  if it has no neighbors.  If nbr_list = NULL, we assume everyone is a neighbor to each other
#' @param inf_inds indices of infectious agents.  indexing starts at 0
#' @param closure_thresh list of length closure_vars where each entry of the list is the closure threshold percentage for the corresponding var
#' @param closure_vars which variables in env_df we will close on
#' @param closure_inds_list list of lists.  first level corresponds to the variable which we are looking at.  second level gives the category of hte variable we are looking at.  in the second level list, there are the indices of all the agents who belong to the same variable category
#' @param closure_time list of length closure_vars where each entry of the list is the time we will close the structure once closure is implemented
#' @param closure_max_T list of length closure_vars where each entry of the list is the max time we will close the structure
#' @return list with new nbr_list, updated do_closure, and updated closure_time, 
update_closure <- function(nbr_list,
                                           inf_inds,
                                           closure_thresh,
                                           closure_vars,
                                           closure_inds_list = NULL,
                                           closure_time = list(0),
                           closure_max_T = list(1)){

    ## Loop over closure vars
    ## Loop over each variable category
    ## Check if closure_time > closure_max_T for variable category
    ## If YES:
    ##    1.  re-add those closure_inds_list back to those indices
    ##    2. turn do_closure to FALSE
    ##
    ## ELSE IF closure_time > 0
    ##   1.  Update closure_time
    ##   2. There is no change, building is still closed.  return same nbr_list
    ##  ELSE closure_time = 0
    ##    1. Update percent infectious for variable category
    ##    IF closure_pct >= closure_thresh
    ##      a. remove closure_inds_list variable category indices from each others nbr_list
    ##      b. update closure_time
    ##    ELSE closure_pct < closure_thresh
    ##       a. Return same nbr_list

    if(length(closure_vars) == 0) stop("closure_vars is length 0")
    for(ii in 1:length(closure_vars)){
        closure_var <- closure_vars[ii]
        n_clos_var_cats <- length(closure_thresh[[ii]])
        for(jj in 1:nclos_var_cats){
            cat_inds <- closure_inds_list[[ii]][[jj]]
            ct <- closure_time[[ii]][jj]
            cT <- closure_max_T[[ii]][jj]
            cThresh <- closure_thresh[[ii]][jj]
            if(ct > cT){
                new_nbr_list <- add_closure_nbr(nbr_list, cat_inds)
                do_closure <- FALSE
            } else if(ct > 0){
                closure_time[[ii]][jj] <- ct + 1
            } else{
                cur_inf_pct <- sum(inf_inds %in% cat_inds) / length(cat_inds)
                if(cur_inf_pct > cThresh){
                    new_nbr_list <- remove_closure_nbr(nbr_list, cat_inds)
                    closure_time[[ii]][jj] <- ct + 1
                }
            }
                
        }

    }
    closure_list <- list(new_nbrs = new_nbr_list,
                         do_closure = do_closure,
                         new_closure_time = closure_time)
    return(closure_list)
    }
