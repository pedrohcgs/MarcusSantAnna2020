#' @title Group-Time Average Treatment Effects for Panel Data and no covariates
#'
#' @description \code{unc_att_gt_panel} computes average treatment effects in DID
#'  setups where there are more than two periods of data and allowing for
#'  treatment to occur at different points in time and allowing for
#'  treatment effect heterogeneity and dynamics. Parallel trends hold without covariates.
#'  See Callaway and Sant'Anna (2019) for a detailed description.
#'
#' @param y_name The name of the outcome variable
#' @param t_name The name of the column containing the time periods
#' @param id_name The individual (cross-sectional unit) id name
#' @param first_treat_name The name of the variable in \code{data} that
#'  contains the first period when a particular observation is treated.
#'  This should be a positive number for all observations in treated groups.
#'  It defines which "group" a unit belongs to.  It should be 0 for units
#'  in the untreated group.
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
unc_att_gt_panel <- function(y_name,
                             t_name,
                             id_name,
                             first_treat_name,
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

  #-----------------------------------------------------------------------------
  # Compute all ATT(g,t)
  #-----------------------------------------------------------------------------
  if(comparison_group == "notyettreated_all"){
    results <- compute_unc_att_gt_panel_ny_all(dp)
  } else {
    results <- compute_unc_att_gt_panel(dp)
  }

  # extract ATT(g,t) and influence functions
  #attgt.list <- results$attgt.list
  #inffunc <- results$inffunc

  # process results
  attgt.results <- process_attgt1(results)
  group <- attgt.results$group
  att <- attgt.results$att
  tt <- attgt.results$tt
  inffunc1 <- attgt.results$inf.func


  # estimate variance
  # this is analogous to cluster robust standard errors that
  # are clustered at the unit level
  n <- dp$n
  V <- t(inffunc1)%*%inffunc1/n
  se <- sqrt(Matrix::diag(V)/n)

  # if clustering along another dimension...we require using the
  # bootstrap (in principle, could come up with analytical standard
  # errors here though)
  if ( (length(cluster_vars) > 0) & !bstrap) {
    warning("clustering the standard errors requires using the bootstrap, resulting standard errors are NOT accounting for clustering")
  }


  # bootstrap variance matrix
  if (bstrap) {
    bout <- did::mboot(inffunc1, DIDparams=dp)
    bres <- bout$bres
    #V <- bout$V
    se <- bout$se
  }


  #-----------------------------------------------------------------------------
  # compute confidence intervals / bands
  #-----------------------------------------------------------------------------

  # critical value from N(0,1), for pointwise
  cval <- stats::qnorm(1-alp/2)

  cband <- TRUE
  # in order to get uniform confidencs bands
  # HAVE to use the bootstrap
  if (bstrap){
    if (cband) {
      # for uniform confidence band
      # compute new critical value
      # see paper for details
      bSigma <- apply(bres, 2,
                      function(b) (stats::quantile(b, .75, type=1, na.rm = T) -
                                     stats::quantile(b, .25, type=1, na.rm = T))/(stats::qnorm(.75) - stats::qnorm(.25)))
      # sup-t confidence band
      bT <- apply(bres, 1, function(b) max( abs(b/bSigma)))
      cval <- stats::quantile(bT, 1-alp, type=1, na.rm = T)
    }
  }

  # Return this list
  return(did::MP(group=group, t=tt, att=att, V_analytical=V, se = se,
                 c=cval, inffunc=inffunc1, n=n, W=NULL, Wpval=NULL, alp = alp, DIDparams=dp))

}
