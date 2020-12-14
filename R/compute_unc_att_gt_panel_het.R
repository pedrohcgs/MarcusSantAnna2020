#' @title Compute Group-Time Average Treatment Effects for Panel Data and no covariates: het case
#'
#' @description \code{compute_unc_att_gt} does the main work for computing
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

#' @export
compute_unc_att_gt_panel_het <- function(dp, het_val){
  #-----------------------------------------------------------------------------
  if(! dp$panel) stop("compute_unc_att_gt only works with Panel Data.")
  #-----------------------------------------------------------------------------
  # unpack DIDparams
  #-----------------------------------------------------------------------------
  data <- dp$data
  yname <- dp$yname
  tname <- dp$tname
  idname <- dp$idname
  printdetails <- dp$print_details
  control.group <- dp$control_group
  first.treat.name <- dp$gname
  n  <- dp$n
  nT <- dp$nT
  nG <- dp$nG
  tlist <- dp$tlist
  glist <- dp$glist


  nevertreated <- (control.group[1] == "nevertreated")
  notyettreated <- (control.group[1] == "notyettreated")
  notyettreated_all <- (control.group[1] == "notyettreated_all")
  #-----------------------------------------------------------------------------
  # will populate with all ATT(g,t)
  attgt.list <- list()

  # place holder in lists
  counter <- 1

  # 3-dimensional array which will store influence function
  # across groups and times
  inffunc <- array(data = 0, dim = c(nG, nT, n))




  # loop over groups
  for (g in 1:nG) {
    ## set an index for the largest pre-treatment period
    pret <- utils::tail(which(tlist < glist[g]),1)
    ## print a error message if there are no pre-treatment period
    if (length(pret) == 0) {
      warning(paste0("There are no pre-treatment periods for the group first treated at ", glist[g]))
      # if there are not pre-treatment periods, code will
      # break, jump out of this loop
      break
    }
    # loop over time periods
    for (t in 1:nT) {
      ## print the details of which iteration we are on
      if (printdetails) {
        cat(paste("current period:", tlist[t]), "\n")
        cat(paste("current group:", glist[g]), "\n")
        cat(paste("Latest pre-treatment period is", tlist[pret]), "\n")
      }

      ## --------------------------------------------------------
      ## results for the case with panel data
      post.treat <- 1 * (glist[g] <= tlist[t])

      if (tlist[pret] == tlist[t]){
        attgt.list[[counter]] <- list(att = 0, group = glist[g], year = tlist[t], post = post.treat)
        inffunc[g,t,] <- 0
      }
      else {
        ## get dataset with current period and latest pre-treatment period
        disdat <- data[((data[, tname] == tlist[t]) | (data[,tname]==tlist[pret])),]
        ## set up control group

        if(nevertreated == TRUE){
          disdat$C <- 1 * (disdat[, first.treat.name] == 0)
        }


        if(notyettreated == TRUE){
          disdat$C <- 1 * ((disdat[, first.treat.name] == 0) +
                             (disdat[, first.treat.name] > max(disdat[, tname]))) *
            (disdat[,first.treat.name] != glist[g])
        }

        ## set up dummy for particular treated group
        disdat$G <- 1 * (disdat[,first.treat.name] == glist[g])

        # transform  disdat it into "cross-sectional" data where one of the columns contains the change in the outcome
        ## over time. dy is computed as latest year - earliest year. We then keep the y of earliest year
        disdat <- suppressWarnings(BMisc::panel2cs(disdat, yname, idname, tname))

        ## drop missing factors
        disdat <- base::droplevels(disdat)

        ## give short names for data in this iteration
        G <- disdat$G
        C <- disdat$C
        dy <- disdat$dy * ((-1)^(1 + post.treat))
        w <- (disdat$w1) * het_val + (disdat$w0) * (1 - het_val)
        #w <- disdat$w

        # The adjustment above in dy is necessary to ensure that the dy has the right sign if post.tread=0
        # since disdat compute Y_last_date - Y_early_date as dy,
        # but with pre-treat it should be Y_early_date - Y_last_date,
        # as last_date is the "pre-treatment period" g-1

        ## set up weights
        attw <- w * G / mean(w * G)
        attw2a <- w * C
        attw2 <- attw2a / mean(attw2a)
        att <- mean((attw - attw2)*dy)

        if(is.na(att)) att <- NA
        if(is.nan(att)) att <- NA

        ## save results for this iteration
        attgt.list[[counter]] <- list(att = att, group = glist[g], year = tlist[t], post = post.treat)

        ## --------------------------------------------
        ## get the influence function

        ## weights
        wg <- w * G/mean(w * G)
        wc1 <- w * C
        wc <- wc1 / mean(wc1)

        ## influence function for treated group
        psig <- wg * (dy - mean(wg * dy))
        # influence function for the control group
        psic <- wc * (dy - mean(wc * dy))

        ## save the influence function as the difference between
        ## the treated and control influence functions;
        ## we save this as a 3-dimensional array
        ## and then process afterwards
        infl.att <- (psig - psic)
        if(base::anyNA(infl.att)) infl.att <- 0
        if(base::any(is.nan(infl.att))) infl.att <- 0

        inffunc[g, t,] <- infl.att
      }

      counter <- counter+ 1
    }
  }


  list(attgt.list = attgt.list, inffunc=inffunc)
}
