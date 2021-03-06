% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_unc_att_gt_panel_ny_all_het.R
\name{compute_unc_att_gt_panel_ny_all_het}
\alias{compute_unc_att_gt_panel_ny_all_het}
\title{Compute ATT(g,t) for Panel Data and no covariates using all available groups}
\usage{
compute_unc_att_gt_panel_ny_all_het(dp, het_val)
}
\arguments{
\item{dp}{A DIDparams object}

\item{het_val}{Either 0 or 1 (value for the subpopulation)}
}
\value{
a list with length equal to the number of groups times the
 number of time periods; each element of the list contains an
 object that contains group-time average treatment effect as well
 as which group it is for and which time period it is for. It also exports
 the influence function which is used externally to compute
 standard errors.
}
\description{
\code{compute_unc_att_gt_panel_ny_all} does the main work for computing
 multiperiod group-time average treatment effects
}
\keyword{internal}
