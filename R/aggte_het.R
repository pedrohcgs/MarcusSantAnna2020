#' @title Differentce of Aggregate ATT's across two subpopulations: Panel data and no covariates
#'
#' @description A function to take group-time average treatment effects
#'  and aggregate them into a smaller number of parameters.  There are
#'  several possible aggregations including "simple", "dynamic", "selective",
#'  and "calendar."
#'
#' @param MP an MP object (i.e., the results of the \code{unc_att_gt_panel_het} method)
#' @param type Which type of aggregated treatment effect parameter to compute.
#'   The default is "simple" (this just computes a weighted average of all
#'   group-time average treatment effects with weights proportional to group
#'   size).  Other options are "dynamic" (this computes average effects across
#'   different lengths of exposure to the treatment and is similar to an
#'   "event study"; here the overall effect averages the effect of the
#'   treatment across all positive lengths of exposure); "selective" (this
#'   computes average treatment effects across different groups; here
#'   the overall effect averages the effect across different groups); and
#'   "calendar" (this computes average treatment effects across different
#'   time periods; here the overall effect averages the effect across each
#'   time period).
#' @param balance.e If set (and if one computes dynamic effects), it balances
#'  the sample with respect to event time.  For example, if \code{balance.e=2},
#'  \code{aggte} will drop groups that are not exposed to treatment for
#'  at least three periods. (the initial period when \code{e=0} as well as the
#'  next two periods when \code{e=1} and the \code{e=2}).  This ensures that
#'  the composition of groups does not change when event time changes.
#'
#' @return AGGTEobj
#' @export
#'
aggte_het <- function(MP, type="simple", balance.e=NULL) {

  MP1 <- MP
  MP1$att <- MP1$att_het1
  MP1$inffunc <- MP1$inffunc1_het1
  MP1$DIDparams$data$w <- MP1$DIDparams$data$w1

  aggte_het1 <- did::compute.aggte(MP1, type, balance.e, na.rm = T)

  MP0 <- MP
  MP0$att <- MP0$att_het0
  MP0$inffunc <- MP0$inffunc1_het0
  MP0$DIDparams$data$w <- MP0$DIDparams$data$w0

  aggte_het0 <- did::compute.aggte(MP0, type, balance.e, na.rm = T)

  inf.function = NULL

  if (type == "dynamic"){
    # Select those with event time within the range of both
    egt_het1 <- aggte_het1$egt
    egt_het0 <- aggte_het0$egt
    egt_all <- base::intersect(egt_het1, egt_het0)
    #e_select <- (egt_all>=mine) & ((egt_all<=maxe))
    #egt_all <- egt_all[e_select]
    sel_het1 <- which(egt_het1 %in% egt_all)
    sel_het0 <- which(egt_het0 %in% egt_all)

    att_select <-  aggte_het1$att.egt[sel_het1] - aggte_het0$att.egt[sel_het0]

    inf.function <- aggte_het1$inf.function$dynamic.inf.func.e[, sel_het1] -
      aggte_het0$inf.function$dynamic.inf.func.e[, sel_het0]

    out <- list(dif = did::AGGTEobj(overall.att = aggte_het1$overall.att -  aggte_het0$overall.att,
                                    overall.se = NULL,
                                    type = type,
                                    egt = egt_all,
                                    att.egt = att_select,
                                    se.egt = NULL,
                                    crit.val.egt = NULL,
                                    inf.function = inf.function),
                aggte_het0 = aggte_het0,
                aggte_het1 = aggte_het1)


  } else {
    out <- list(aggte_het0 = aggte_het0,
                aggte_het1 = aggte_het1)

  }
  return(out)

}
