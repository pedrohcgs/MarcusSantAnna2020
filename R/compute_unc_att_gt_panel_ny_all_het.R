#' @title Compute ATT(g,t) for Panel Data and no covariates using all available groups
#'
#' @description \code{compute_unc_att_gt_panel_ny_all} does the main work for computing
#'  multiperiod group-time average treatment effects
#'
#' @param dp A DIDparams object
#' @param het_val Either 0 or 1 (value for the subpopulation)
#'
#' @return a list with length equal to the number of groups times the
#'  number of time periods; each element of the list contains an
#'  object that contains group-time average treatment effect as well
#'  as which group it is for and which time period it is for. It also exports
#'  the influence function which is used externally to compute
#'  standard errors.
#'
#' @keywords internal
#'
#' @export
compute_unc_att_gt_panel_ny_all_het <- function(dp, het_val) {

  #-----------------------------------------------------------------------------
  # unpack DIDparams
  #-----------------------------------------------------------------------------
  data <- dp$data
  yname <- dp$yname
  tname <- dp$tname
  idname <- dp$idname
  printdetails <- dp$printdetails
  control.group <- dp$control.group
  first.treat.name <- dp$first.treat.name
  n  <- dp$n
  nT <- dp$nT
  nG <- dp$nG
  tlist <- dp$tlist
  glist <- dp$glist

  #nevertreated <- (control.group[1] == "nevertreated")
  #notyettreated <- (control.group[1] == "notyettreated")
  notyettreated_all <- (control.group[1] == "notyettreated_all")
  #-----------------------------------------------------------------------------
  if(!notyettreated_all) stop("This code is only suitable for comparison_group = not_yet_all.")

  #-----------------------------------------------------------------------------
  # main computations
  #-----------------------------------------------------------------------------

  # will populate with all att(g,t)
  attgt.list <- list()

  # place holder in lists
  counter <- 1

  # number of time periods
  tlist.length <- length(tlist)

  # 3-dimensional array which will store influence function
  # across groups and times
  inffunc <- array(data=0, dim=c(nG,nT,n))

  # loop over groups
  for (g in 1:nG) {
    # These are useful for computing the recursive components of the ATT(g,t)
    att_previous = 0
    infl.att_previous = 0

    # loop over time periods
    for (t in 1:(tlist.length-1)) {


      # set pre-treatment time period (this is updated later
      # if g <= t (i.e. for "already treated" groups)
      pret <- t


      # code to update pre-treatment time periods
      if (glist[g]<=tlist[(t + 1)]) {

        # set an index for the pretreatment period
        # this recovers the right pre-treatment period for this group
        # it is the most recent pre-treatment period (g-1)
        pret <- utils::tail(which(tlist < glist[g]),1)

        # stop if there are no pre-treatment period
        if (length(pret) == 0) {
          warning(paste0("There are no pre-treatment periods for the group first treated at ", glist[g]))
          break
        }

      }

      # print the details of which iteration we are on
      if (printdetails) {
        cat(paste("current period:", tlist[(t + 1)]), "\n")
        cat(paste("current group:", glist[g]), "\n")
        cat(paste("set pre-treatment period to be", tlist[pret]), "\n")
      }

      #-----------------------------------------------------------------------------
      # results for the case with panel data
      #-----------------------------------------------------------------------------
      # post treatment dummy variable
      post.treat <- 1 * (glist[g] <= tlist[(t + 1)])
      # Get instantaneous treatment indicator
      #inst.treat <- 1 * (glist[g] == tlist[t + 1])
      #-----------------------------------------------------------------------------
      if (notyettreated_all) pret <- t

      # get dataset with current period and pre-treatment period
      # Pre-treatment period here is t-1
      disdat <- data[(data[,tname]==tlist[t+1] | data[,tname]==tlist[pret]),]

      # set up control group
      # use "not yet treated as control"
      # that is, never treated + units that are eventually treated,
      # but not treated by the current period (excluding the group g that is being uses!)
      disdat$C <- 1 * ((disdat[,first.treat.name] == 0) |
                         ((disdat[,first.treat.name] > tlist[t+1]) &
                            (disdat[,first.treat.name] != glist[g])))

      # set up dummy for particular treated group
      disdat$G <- 1*(disdat[,first.treat.name] == glist[g])


      # transform  disdat it into "cross-sectional" data where one of the columns
      # contains the change in the outcome over time.
      # dy is computed as latest year - earliest year. "Y" is outcome
      # in the pre period, "yt1" is outcome in the post period
      disdat <- BMisc::panel2cs(disdat, yname, idname, tname)

      # drop missing factors
      disdat <- base::droplevels(disdat)

      # give short names for data in this iteration
      G <- disdat$G
      C <- disdat$C
      dy <- disdat$dy
      w <- (disdat$w1) * het_val + (disdat$w0) * (1 - het_val)
      #w <- disdat$w

      #-----------------------------------------------------------------------------
      # code for actually computing ATT(g,t)
      #-----------------------------------------------------------------------------
      ## set up weights and compute E[delta_Y|Gg=1] - E[delta_Y|Gg=1]
      att_treated_w <- w * G / mean(w * G)
      att_comp_w <- w * C / mean(w * C)
      #-----------------------------------------------------------------------------
      att_treated <- mean(att_treated_w * dy)
      att_comp <- mean(att_comp_w * dy)

      # get the ATT(g,t), which is now recursive (sum from s=g to t)
      att <- att_treated - att_comp + (att_previous)

      if(is.na(att)) att <- NA
      if(is.nan(att)) att <- NA

      ## save results for this iteration
      attgt.list[[counter]] <- list(att = att, group = glist[g], year = tlist[t+1], post = post.treat)

      ## --------------------------------------------
      ## get the influence function
      ## influence function for treated group
      psig <- att_treated_w * (dy - att_treated)
      # influence function for the control group
      psic <- att_comp_w * (dy - att_comp)

      ## save the influence function as the difference between
      ## the treated and control influence functions;
      ## we save this as a 3-dimensional array
      ## and then process afterwards
      infl.att <- psig - psic + (infl.att_previous)
      if(base::anyNA(infl.att)) infl.att <- 0
      if(base::any(is.nan(infl.att))) infl.att <- 0

      inffunc[g, (t + 1),] <- infl.att
      ## --------------------------------------------
      # Save these info that will be used for t>g
      # Save these info that will be used for t>g
      att_previous <- att * post.treat
      infl.att_previous <- infl.att * post.treat
      ## --------------------------------------------
      # update counter
      counter <- counter + 1
    } # end looping over t
  } # end looping over g

  return(list(attgt.list = attgt.list, inffunc = inffunc))
}
