


process_attgt1 <- function(attgt.results.list) {
  attgt.list <- attgt.results.list$attgt
  inffunc <- attgt.results.list$inffunc
  nG <- length(unique(unlist(BMisc::getListElement(attgt.list, "group"))))
  nT <- length(unique(unlist(BMisc::getListElement(attgt.list, "year"))))
  # pick up number of observations from the influence function

  n <- length(inffunc[1,1,])
  # create vectors to hold the results
  group <- c()
  att <- c()
  tt <- c()
  i <- 1


  inffunc1 <- matrix(0, ncol=nG*nT, nrow=n)

  # populate result vectors and matrices
  for (f in 1:nG) {
    for (s in 1:nT) {
      group[i] <- attgt.list[[i]]$group
      tt[i] <- attgt.list[[i]]$year
      att[i] <- attgt.list[[i]]$att
      inffunc1[,i] <- inffunc[f,s,]
      i <- i+1
    }
  }

  list(group=group, att=att, tt=tt, inf.func=inffunc1)
}
