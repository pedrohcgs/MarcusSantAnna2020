#' @title Difference of ATT(g,t)'s across two subpopulations: Panel data and no covariates
#'
#' @description \code{unc_att_gt_panel_het2} computes the difference of average treatment effects across two subpopulations
#' in DID setups where there are more than two periods of data and
#' allowing for treatment to occur at different points in time. Here, we assume same trends
#' between the two supopulations. See Marcus ans Sant'Anna (2020) for a detailed description.
#'
#'
#' @param y_name The name of the outcome variable
#' @param t_name The name of the column containing the time periods
#' @param id_name The individual (cross-sectional unit) id name
#' @param first_treat_name The name of the variable in \code{data} that
#'  contains the first period when a particular observation is treated.
#'  This should be a positive number for all observations in treated groups.
#'  It defines which "group" a unit belongs to.  It should be 0 for units
#'  in the untreated group.
#' @param het_name The name of the column containing the binary categories for heterogeneity
#' @param weights_name The name of the column containing the sampling weights.
#'  If not set, all observations have same weight.
#' @param cluster_vars A vector of variables to cluster on.  At most, there
#'  can be two variables (otherwise will throw an error) and one of these
#'  must be the same as idname which allows for clustering at the individual
#'  level.
#' @param comparison_group Which units to use the comparison group.
#'  The default is \code{group="never"} which sets the comparison group
#'  to be the group of units that never participate in the
#'  treatment.  This group does not change across groups or
#'  time periods. The other option is to set
#'  \code{group="not_yet"}.  In this case, the comparison group
#'  is set to the group of units that have not yet participated
#'  in the treatment in that time period.  This includes all
#'  never treated units, but it includes additional units that
#'  eventually participate in the treatment, but have not
#'  participated yet.
#'
#' @param data The name of the data.frame that contains the data
#' @param alp the significance level, default is 0.05
#' @param bstrap Boolean for whether or not to compute standard errors using
#'  the multiplier boostrap.  If standard errors are clustered, then one
#'  must set \code{bstrap=TRUE}. Default is \code{TRUE} (in addition, cband
#'  is also by default \code{TRUE} indicating that uniform confidence bands
#'  will be returned.  If bstrap is \code{FALSE}, then analytical
#'  standard errors are reported.
#' @param nboot The number of boostrap iterations to use.  The default is 999,
#'  and this is only applicable if \code{bstrap=TRUE}.
#' @param print_details Boolean for showing detailed results or not
#' @param parallel Boolean for whether or not to use parallel processing
#'  (not implemented yet)
#' @param cores The number of cores to use for parallel processing
#'  (not implemented yet)
#' @references Callaway, Brantly and Sant'Anna, Pedro H. C.. "Difference-in-Differences with Multiple Time Periods and an Application on the Minimum Wage and Employment." Working Paper <https://ssrn.com/abstract=3148250> (2019).
#' @references Marcus, Michelle and Sant'Anna, Pedro H. C.. "The Role of Parallel Trends in Event Study Settings: An Application to Environmental Economics." Working Paper <https://to_be_added.com> (2020).
#'
#'
#' @return an \code{MP} object containing all the results for group-time average
#'  treatment effects
#'
#' @export
#'
unc_att_gt_panel_het2 <- function(y_name,
                                 t_name,
                                 id_name,
                                 first_treat_name,
                                 het_name,
                                 weights_name = NULL,
                                 cluster_vars = NULL,
                                 comparison_group = c("never", "not_yet", "not_yet_all"),
                                 data,
                                 alp = 0.05,
                                 bstrap = TRUE,
                                 nboot = 999,
                                 print_details = TRUE,
                                 parallel = FALSE,
                                 cores = 1) {

  # Map some parameters
  if(comparison_group=="never") comparison_group <- "nevertreated"
  if(comparison_group=="not_yet") comparison_group <- "notyettreated"
  if(comparison_group=="not_yet_all") comparison_group <- "notyettreated_all"


  # het if null
  if(is.null(het_name)) {
    stop("Please specifiy 'het_name'. If het_name=NULL, use 'att_gt' instead of 'att_gt_het'.")
  }

  if(!is.character(het_name)) stop(" 'het_name' must be a character")

  het <- as.vector(data[, het_name])

  #Check for missing values in het_name (to avoid needing to recode things up)
  if (!all(!is.na(het))) stop(" The column of 'het_name' can not contain missing values.")

  het.dim <- base::nrow(unique(het))
  if( het.dim!=2 )  {
    stop("'het_name' must be a binary variable.")
  }

  if(sum((het == 1) + (het == 0)) != base::nrow(het))  {
    stop("'het_name' must be a binary (0-1) variable.")
  }


  # this is a DIDparams object
  dp <- did::pre_process_did(yname = y_name,
                             tname = t_name,
                             idname = id_name,
                             gname = first_treat_name,
                             xformla = NULL,
                             data = data,
                             panel = TRUE,
                             allow_unbalanced_panel = FALSE,
                             control_group = comparison_group,
                             weightsname = weights_name,
                             alp = alp,
                             bstrap = bstrap,
                             cband = TRUE,
                             biters = nboot,
                             clustervars = cluster_vars,
                             est_method = "ipw",
                             print_details = print_details,
                             pl = parallel,
                             cores = cores
  )

  het <- dp$data[, het_name]

  dp$data$w1 <- dp$data$w * (het == 1)
  dp$data$w0 <- dp$data$w * (het == 0)

  #-----------------------------------------------------------------------------
  # Compute all ATT(g,t)
  #-----------------------------------------------------------------------------
  if(comparison_group == "notyettreated_all"){
    results_het1 <- compute_unc_att_gt_panel_ny_all_het2(dp, het_val = 1)
    results_het0 <- compute_unc_att_gt_panel_ny_all_het2(dp, het_val = 0)
  } else {
    results_het1 <- compute_unc_att_gt_panel_het2(dp, het_val = 1)
    results_het0 <- compute_unc_att_gt_panel_het2(dp, het_val = 0)
  }
  #-----------------------------------------------------------------------------
  # Extract info for het ==1
  # extract ATT(g,t) and influence functions
  #attgt.list_het1 <- results_het1$attgt.list
  #inffunc_het1 <- results_het1$inffunc
  # process results
  attgt.results_het1 <- process_attgt1(results_het1)
  att_het1 <- attgt.results_het1$att
  inffunc1_het1 <- attgt.results_het1$inf.func

  group_het1 <- attgt.results_het1$group
  tt_het1 <- attgt.results_het1$tt
  #-----------------------------------------------------------------------------
  # Extract info for het ==0
  # extract ATT(g,t) and influence functions
  #attgt.list_het0 <- results_het0$attgt.list
  #inffunc_het0 <- results_het0$inffunc
  # process results
  attgt.results_het0 <- process_attgt1(results_het0)
  att_het0 <- attgt.results_het0$att
  inffunc1_het0 <- attgt.results_het0$inf.func

  group_het0 <- attgt.results_het0$group
  tt_het0 <- attgt.results_het0$tt
  #-----------------------------------------------------------------------------
  if(!all.equal(group_het0, group_het1)) stop("Something is wrong here as the groups should coincide")
  if(!all.equal(tt_het0, tt_het1)) stop("Something is wrong here as the times should coincide")
  #-----------------------------------------------------------------------------
  out <- list(
    group = group_het0,
    t = tt_het0,
    att_het1 = att_het1,
    att_het0 = att_het0,
    inffunc1_het1 = inffunc1_het1,
    inffunc1_het0 = inffunc1_het0,
    alp = alp,
    DIDparams = dp)

  class(out) <- "MP"



  # Return this list
  return(out)

}
