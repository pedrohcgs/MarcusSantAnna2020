# THIS IS FROM BRANT
#' @title Compute Aggregated Treatment Effect Parameters
#'
#' @description Does the heavy lifting on computing aggregated group-time
#'  average treatment effects
#'
#' @inheritParams did::att_gt
#' @inheritParams did::aggte
#'
#' @return \code{AGGTEobj} object - check did package
#'
#'
#' @export
compute_aggte_new <- function(MP, type="simple", balance.e=NULL) {

  #-----------------------------------------------------------------------------
  # unpack MP object
  #-----------------------------------------------------------------------------
  # load parameters
  group <- MP$group
  t <- MP$t
  att <- MP$att
  inffunc1 <- MP$inffunc
  n <- MP$n

  dp <- MP$DIDparams
  first.treat.name <- dp$first.treat.name
  clustervars <- dp$clustervars
  data <- dp$data
  tname <- dp$tname
  idname <- dp$idname
  bstrap <- dp$bstrap
  biters <- dp$biters
  alp <- dp$alp
  cband <- dp$cband
  tlist <- dp$tlist
  glist <- dp$glist
  panel <- dp$panel

  # data from first period
  ifelse(panel,
         dta <- data[ data[,tname]==tlist[1], ],
         dta <- data
  )
  #-----------------------------------------------------------------------------
  # data organization and recoding
  #-----------------------------------------------------------------------------

  # do some recoding to make sure time periods are 1 unit apart
  # and then put these back together at the end
  originalt <- t
  originalgroup <- group
  originalglist <- glist
  originaltlist <- tlist
  uniquet <- seq(1,length(unique(t)))
  # function to switch from "new" t values to  original t values
  t2orig <- function(t) {
    unique(c(originalt,0))[which(c(uniquet,0)==t)]
  }
  # function to switch between "original"
  #  t values and new t values
  orig2t <- function(orig) {
    c(uniquet,0)[which(unique(c(originalt,0))==orig)]
  }
  t <- sapply(originalt, orig2t)
  group <- sapply(originalgroup, orig2t)
  glist <- sapply(originalglist, orig2t)

  tlist <- unique(t)
  maxT <- max(t)

  # Set the weights
  weights.ind  <-  dta$w

  ## some variables used throughout
  # Ever treated only among the units we actually compute the ATT(g,t)
  ever.treated <- 1 * (dta[, first.treat.name]>0)
  mean.w.ever.treated <- mean(weights.ind * ever.treated)

  # Probability of being in group g (among the relevant ever-treated groups)
  pg <- sapply(originalglist,
               function(g) {
                 mean(weights.ind * ever.treated * (dta[, first.treat.name] == g))/
                   mean(weights.ind * ever.treated)
               })

  # length of this is equal to number of groups
  pgg <- pg

  # same but length is equal to the number of ATT(g,t)
  pg <- pg[match(group, glist)]

  # which group time average treatment effects are post-treatment
  keepers <- which(group <= t)

  # n x 1 vector of group variable
  G <-  unlist(lapply(dta[,first.treat.name], orig2t))

  #-----------------------------------------------------------------------------
  # Compute the simple ATT summary
  #-----------------------------------------------------------------------------

  if (type == "simple") {

    # simple att
    # averages all post-treatment ATT(g,t) with weights
    # given by group size
    simple.att <- sum(att[keepers] * pg[keepers])/(sum(pg[keepers]))
    #-----------------------------------------------------------------------------
    #Now get influence functions

    # Estimation effect coming from P(G=g| Ever treated)
    # Part 1: est effect from P(G=g, ever treated = 1) treating P(ever treated) as known
    simple.oif1 <- sapply(keepers,
                          function(k) {
                            (weights.ind *(G == group[k]) - mean( weights.ind * (G == group[k]))) /
                              mean.w.ever.treated
                          })

    # Part 2: est effect from  P(ever treated) treating P(G=g, ever treated = 1) as known
    simple.oif2 <- sapply(keepers,
                          function(j) {
                            (mean(weights.ind * (G == group[j])) / (mean.w.ever.treated^2)) *
                              (weights.ind * ever.treated - mean.w.ever.treated)
                          })

    # Estimation effect from numerator
    simple.oif <- (simple.oif1 - simple.oif2)/(sum(pg[keepers]))
    #Estimation effect from denominator of the weights (normalization)
    simple.oif3 <- base::rowSums(simple.oif) %*%  t(matrix(pg[keepers]/sum(pg[keepers])))
    #Estimation effect from estimated weights (in total)
    simple.wif <- simple.oif - simple.oif3

    # get the overall influence function
    simple.if <- get_agg_inf_func(att = att,
                                  inffunc1 = inffunc1,
                                  whichones = keepers,
                                  weights.agg = pg[keepers]/sum(pg[keepers]),
                                  wif = simple.wif)

    # get standard errors from overall influence function
    simple.se <- getSE(simple.if, dp)

    return(did::AGGTEobj(overall.att = simple.att,
                         overall.se = simple.se,
                         type = type,
                         inf.function = list(simple.att = simple.if)))
  }

   #-----------------------------------------------------------------------------
  # Compute the event-study estimators
  #-----------------------------------------------------------------------------

  if (type == "dynamic") {

    # event times
    # this looks at all available event times
    # note: event times can be negative here.
    # note: event time = 0 corresponds to "on impact"
    #eseq <- unique(t-group)
    eseq <- unique(originalt - originalgroup)
    #eseq <- eseq[order(eseq)]
    maxe <- max(eseq)
    eseq <- 0 : maxe

    # if the user specifies balance.e, then we are going to
    # drop some event times and some groups; if not, we just
    # keep everything (that is what this variable is for)
    include.balanced.gt <- rep(TRUE, length(originalgroup))

    # if we balance the sample with respect to event time
    if (!is.null(balance.e)) {
      eseq <- eseq[ (eseq <= balance.e) & (eseq >= balance.e - t2orig(maxT) + t2orig(1))]
      include.balanced.gt <- (t2orig(maxT) - originalgroup >= balance.e)
    }

    # these are not currently used, but if we want to trim
    # out some lengths of exposure, we could use this
    # eseq <- eseq[ (eseq >= mine) & (eseq <= maxe) ]
    # note that they would still be included in estimating overall effects

    # compute atts that are specific to each event time
    dynamic.att.e <- sapply(eseq, function(e) {
      # keep att(g,t) for the right g&t as well as ones that
      # are not trimmed out from balancing the sample
      whiche <- which( (originalt - originalgroup == e) & (include.balanced.gt) )
      atte <- att[whiche]
      pge <- pg[whiche]/(sum(pg[whiche]))
      sum(atte*pge)
    })

    # Amount of time each unit is treated
    time.treated <- ((max(originalt)) - dta[,first.treat.name] + 1) * (dta[, first.treat.name] > 0 )

    # Get the influence function of the event-study estimators (post-treatment)
    dynamic.e.inf.f <- sapply(eseq, function(e) {
      whiche <- which((t - group) == e)
      pge <- pg[whiche]/sum(pg[whiche])

      ## some variables used
      atleast.e.treated <- 1 * ((time.treated - 1) >= e)
      mean.w.atleast.e.treated <- mean(weights.ind * atleast.e.treated)

      # Estimation effect coming from P(G=g| treated for at least e periods)
      # Part 1: est effect from P(G=g, treated at least e periods) treating P(treated for at least e periods) as known
      dynamic.oif1 <- sapply(whiche,
                             function(k) ( (weights.ind *(G==group[k]) - mean(weights.ind * (G==group[k]))) /
                                             mean.w.atleast.e.treated
                             )
      )

      # Part 2: est effect from  P(treated for at least e periods) treating P(G=g) as known
      dynamic.oif2 <- sapply(whiche,
                             function(j) ((mean(weights.ind * (G==group[j]))/(mean.w.atleast.e.treated^2)) *
                                            (weights.ind * atleast.e.treated - mean.w.atleast.e.treated)
                             )
      )

      dynamic.oif <- (dynamic.oif1 - dynamic.oif2)

      inffunc1[,whiche] %*% as.matrix(pge) + dynamic.oif %*% as.matrix(att[whiche])

      #getSE(whiche, pge, dynamic.oif )

    })

        # Pre-treatment periods
    pre.treat <- which(group > t)

    mine <- min(t-group)
    eseq.pre <- seq(mine, -1)


    dynamic.att.pre.e <- sapply(eseq.pre, function(e) {
      whiche <- which(t - group == e)
      atte <- att[whiche]
      pge <- pg[whiche]/(sum(pg[whiche]))
      sum(atte*pge)
    })


    n.pre.treat <- (time.treated - max(t) + 1) * (dta[, first.treat.name] > 0)


    dynamic.pre.e.inf.f <- sapply(eseq.pre, function(e) {
      whiche <- which(t - group == e)
      pge <- pg[whiche]/sum(pg[whiche])

      ## some variables used
      atleast.e.pre.treated <- 1 * ((n.pre.treat - 1) <= e) * (dta[,first.treat.name]>0)
      mean.w.atleast.e.pre.treated <- mean(weights.ind * atleast.e.pre.treated)

      # Estimation effect coming from P(G=g| treated for at least e periods)
      # Part 1: est effect from P(G=g) treating P(treated for at least e periods) as known
      dynamic.oif1 <- sapply(whiche,
                             function(k) ( (weights.ind *(G==group[k]) - mean(weights.ind * (G==group[k]))) /
                                             mean.w.atleast.e.pre.treated
                             )
      )

      # Part 2: est effect from  P(treated for at least e periods) treating P(G=g) as known
      dynamic.oif2 <- sapply(whiche,
                             function(j) ((mean(weights.ind * (G==group[j]))/(mean.w.atleast.e.pre.treated^2)) *
                                            (weights.ind * atleast.e.pre.treated - mean.w.atleast.e.pre.treated)
                             )
      )
      dynamic.oif <- (dynamic.oif1 - dynamic.oif2)

      inffunc1[,whiche] %*% as.matrix(pge) + dynamic.oif %*% as.matrix(att[whiche])
    })

    dynamic.inf.func.e <- cbind(dynamic.pre.e.inf.f, dynamic.e.inf.f)
    colnames(dynamic.inf.func.e) <- c(eseq.pre, eseq)

    dynamic.se.e <- sqrt(base::colMeans(dynamic.inf.func.e^2)/n)

    dynamic.crit.val <- NULL


    # Average over all dynamic att(e), e>=0
    dynamic.att <- mean(dynamic.att.e)
    # Influence function of the average over all dynamic att(e), e>=0
    dynamic.if <- rowMeans(dynamic.e.inf.f)



    # get overall average treatment effect
    # by averaging over positive dynamics
    #epos <- (eseq >= 0)
    #dynamic.att <- mean(dynamic.att.e[epos])
    #dynamic.inf.func <- get_agg_inf_func(att = dynamic.att.e[epos],
    #                                     inffunc1 = as.matrix(dynamic.inf.func.e[,epos]),
    #                                     whichones = (1:sum(epos)),
    #                                     weights.agg = (rep(1/sum(epos), sum(epos))),
    #                                     wif = NULL)
    dynamic.se <- getSE(dynamic.if, dp)

    #print(dynamic.if - dynamic.inf.func)

    return(did::AGGTEobj(overall.att = dynamic.att,
                         overall.se = dynamic.se,
                         type = type,
                         egt = eseq,
                         att.egt = c(dynamic.att.pre.e, dynamic.att.e),
                         se.egt = dynamic.se.e,
                         crit.val.egt = dynamic.crit.val,
                         inf.function = list(dynamic.inf.func.e = dynamic.inf.func.e,
                                             dynamic.inf.func = dynamic.if
                                             )))
  }

  #-----------------------------------------------------------------------------
  # Compute the selective treatment timing estimators
  #-----------------------------------------------------------------------------
  if (type == "selective") {

    # get group specific ATTs
    # note: there are no estimated weights here
    selective.att.g <- sapply(glist, function(g) {
      # look at post-treatment periods for group g
      whichg <- which( (group == g) & (g <= t))
      attg <- att[whichg]
      mean(attg)
    })

    # get standard errors for each group specific ATT
    selective.se.inner <- lapply(glist, function(g) {
      whichg <- which( (group == g) & (g <= t))
      inf.func.g <- get_agg_inf_func(att=att,
                                     inffunc1=inffunc1,
                                     whichones=whichg,
                                     weights.agg=pg[whichg]/sum(pg[whichg]),
                                     wif=NULL)
      se.g <- getSE(inf.func.g, dp)
      list(inf.func=inf.func.g, se=se.g)
    })

    # recover standard errors separately by group
    selective.se.g <- unlist(getListElement(selective.se.inner, "se"))

    # recover influence function separately by group
    selective.inf.func.g <- as.matrix(simplify2array(getListElement(selective.se.inner, "inf.func"))[,1,])

    # use multiplier boostrap (across groups) to get critical value
    # for constructing uniform confidence bands
    selective.crit.val <- did::mboot(selective.inf.func.g, dp)$crit.val

    # get overall att under selective treatment timing
    # (here use pgg instead of pg because we can just look at each group)
    selective.att <- sum(selective.att.g * pgg)/sum(pgg)

    # account for having to estimate pgg in the influence function
    selective.wif <- wif(keepers=1:length(glist),
                         pg=pgg,
                         weights.ind=weights.ind,
                         G=G,
                         group=group)

    # get overall influence function
    selective.inf.func <- get_agg_inf_func(att=selective.att.g,
                                           inffunc1=selective.inf.func.g,
                                           whichones=(1:length(glist)),
                                           weights.agg=pgg/sum(pgg),
                                           wif=selective.wif)

    # get overall standard error
    selective.se <- getSE(selective.inf.func, dp)

    return(did::AGGTEobj(overall.att = selective.att,
                         overall.se = selective.se,
                         type = type,
                         egt = originalglist,
                         att.egt = selective.att.g,
                         se.egt = selective.se.g,
                         crit.val.egt = selective.crit.val))

  }




  #-----------------------------------------------------------------------------
  # calendar time effects
  #-----------------------------------------------------------------------------

  if (type == "calendar") {

    # drop time periods where no one is treated yet
    # (can't get treatment effects in those periods)
    minG <- min(group)
    calendar.tlist <- tlist[tlist>=minG]

    # calendar time specific atts
    calendar.att.t <- sapply(calendar.tlist, function(t1) {
      # look at post-treatment periods for group g
      whicht <- which( (t == t1) & (group <= t))
      attt <- att[whicht]
      mean(attt)
    })

    # get standard errors and influence functions
    # for each time specific att
    calendar.se.inner <- lapply(calendar.tlist, function(t1) {
      whicht <- which( (t == t1) & (group <= t))
      wif.t <- wif(keepers=whicht,
                   pg=pg,
                   weights.ind=weights.ind,
                   G=G,
                   group=group)
      inf.func.t <- get_agg_inf_func(att=att,
                                     inffunc1=inffunc1,
                                     whichones=whicht,
                                     weights.agg=pg[whicht]/sum(pg[whicht]),
                                     wif=wif.t)
      se.t <- getSE(inf.func.t, dp)
      list(inf.func=inf.func.t, se=se.t)
    })

    # recover standard errors separately by time
    calendar.se.t <- unlist(getListElement(calendar.se.inner, "se"))

    # recover influence function separately by time
    calendar.inf.func.t <- as.matrix(simplify2array(getListElement(calendar.se.inner, "inf.func"))[,1,])

    # use multiplier boostrap (across groups) to get critical value
    # for constructing uniform confidence bands
    calendar.crit.val <- did::mboot(calendar.inf.func.t, dp)$crit.val

    # get overall att under calendar time effects
    # this is just average over all time periods
    calendar.att <- mean(calendar.att.t)

    # get overall influence function
    calendar.inf.func <- get_agg_inf_func(att=calendar.att.t,
                                          inffunc1=calendar.inf.func.t,
                                          whichones=(1:length(calendar.tlist)),
                                          weights.agg=rep(1/length(calendar.tlist), length(calendar.tlist)),
                                          wif=NULL)

    # get overall standard error
    calendar.se <- getSE(calendar.inf.func, dp)

    return(did::AGGTEobj(overall.att=calendar.att,
                         overall.se=calendar.se,
                         type=type,
                         egt=sapply(calendar.tlist,t2orig),
                         att.egt=calendar.att.t,
                         se.egt=calendar.se.t,
                         crit.val.egt=calendar.crit.val))

  }


}

#-----------------------------------------------------------------------------
# Internal functions for getteing standard errors
#-----------------------------------------------------------------------------

#' @title Compute extra term in influence function due to estimating weights
#'
#' @description A function to compute the extra term that shows up in the
#'  influence function for aggregated treatment effect parameters
#'  due to estimating the weights
#'
#' @param keepers a vector of indices for which group-time average
#'  treatment effects are used to compute a particular aggregated parameter
#' @param pg a vector with same length as total number of group-time average
#'  treatment effects that contains the probability of being in particular group
#' @param weights.ind additional sampling weights (nx1)
#' @param G vector containing which group a unit belongs to (nx1)
#' @param group vector of groups
#'
#' @return nxk influence function matrix
#'
#' @keywords internal
wif <- function(keepers, pg, weights.ind, G, group) {
  # note: weights are all of the form P(G=g|cond)/sum_cond(P(G=g|cond))
  # this is equal to P(G=g)/sum_cond(P(G=g)) which simplifies things here

  # effect of estimating weights in the numerator
  if1 <- sapply(keepers, function(k) {
    (weights.ind * 1*(G==group[k]) - pg[k]) /
      sum(pg[keepers])
  })
  # effect of estimating weights in the denominator
  if2 <- rowSums( sapply( keepers, function(k) {
    weights.ind*1*(G==group[k]) - pg[k]
  })) %*%
    t(pg[keepers]/(sum(pg[keepers])^2))

  # return the influence function for the weights
  if1 - if2
}


#' @title Get an influence function for particular aggregate parameters
#'
#' @title This is a generic internal function for combining influence
#'  functions across ATT(g,t)'s to return an influence function for
#'  various aggregated treatment effect parameters.
#'
#' @param att vector of group-time average treatment effects
#' @param inffunc1 influence function for all group-time average treatment effects
#'  (matrix)
#' @param whichones which elements of att will be used to compute the aggregated
#'  treatment effect parameter
#' @param weights.agg the weights to apply to each element of att[whichones];
#'  should have the same dimension as att[whichones]
#' @param wif extra influence function term coming from estimating the weights;
#'  should be n x k matrix where k is dimension of whichones
#'
#' @return nx1 influence function
#'
#' @keywords internal
get_agg_inf_func <- function(att, inffunc1, whichones, weights.agg, wif=NULL) {
  # enforce weights are in matrix form
  weights.agg <- as.matrix(weights.agg)

  # multiplies influence function times weights and sums to get vector of weighted IF (of length n)
  thisinffunc <- inffunc1[,whichones]%*%weights.agg

  # Incorporate influence function of the weights
  if (!is.null(wif)) {
    thisinffunc <- thisinffunc + wif%*%as.matrix(att[whichones])
  }

  # return influence function
  return(thisinffunc)
}


#' @title Take influence function and return standard errors
#'
#' @description Function to take an nx1 influence function and return
#'  a standard error
#'
#' @param thisinffunc An influence function
#' @inheritParams compute_aggte_new
#'
#' @return scalar standard error
#'
#' @keywords internal
getSE <- function(thisinffunc, DIDparams=NULL) {
  alp <- .05
  bstrap <- FALSE
  if (!is.null(DIDparams)) {
    bstrap <- DIDparams$bstrap
    alp <- DIDparams$alp
    cband <- DIDparams$cband
    n <- length(thisinffunc)
  }

  if (bstrap) {
    bout <- did::mboot(thisinffunc, DIDparams)
    return(bout$se)
  } else {
    return(sqrt( mean((thisinffunc)^2)/n ))
  }
}

