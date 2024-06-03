simu_Report <- function(envProj){
  dropState <- get("dropState", pos = envProj)
  DropMode <- get("DropMode", pos = envProj)
  FreqMode <- get("FreqMode", pos = envProj)
  FixedLR <- get("FixedLR", pos = envProj)
  PrD <- get("PrD", pos = envProj)
  min.freq <- get("min.freq", pos = envProj)
  freqFp <- get("freqFp", pos = envProj)
  mutationFp <- get("mutationFp", pos = envProj)
  gtFp <- get("gtFp", pos = envProj)
  report_locus_set_name <- get("locus_set_name", pos = envProj)

  softVer <- get("softVer", pos = envProj)
  locus_set <- get("locus_set", pos = envProj)
  linkage.loci <- get("linkage.loci", pos = envProj)
  RR <- get("RR", pos = envProj)
  n.simu <- get("n.simu", pos = envProj)

  return_information.simu <- list()
  return_information.simu <- get("return_information.simu", pos = envProj)

  LR.simu <- return_information.simu[[1]]
  calc.time.simu <- return_information.simu[[2]]
  today.simu <- return_information.simu[[3]]
  df <- list()
  df[[1]] <- get("df_Hypo1", pos = envProj)
  df[[2]] <- get("df_Hypo2", pos = envProj)

  simu_items <- c("","min","5th percentile","25th percentile","50th percentile","75th percentile","95th percentile","max")
  report_today.simu <- format(today.simu, "%Y/%b/%d %X")
  report_calc.time.simu <- paste(signif(calc.time.simu[3], digits=3), "sec", collapse="")

  LRprod <- list()
  for (k in 1:2){
    if (!is.na(RR[1])){
      LRprod[[k]] <- sort(apply(LR.simu[[k]][[1]],1,prod) * apply(LR.simu[[k]][[2]],1,prod))
    }
    else{
      LRprod[[k]] <- sort(apply(LR.simu[[k]],1,prod))
    }
  }

  report_simu<- vector("list", length = 2)
  for (k in 1:2){
    report_simu[[k]][[1]] <- paste(simu_items[2], quantile(LRprod[[k]])["0%"], sep=":,")
    report_simu[[k]][[2]] <- paste(simu_items[3], LRprod[[k]][floor(length(LRprod[[k]]) / 20) + 1], sep=":,")
    report_simu[[k]][[3]] <- paste(simu_items[4], quantile(LRprod[[k]])["25%"], sep=":,")
    report_simu[[k]][[4]] <- paste(simu_items[5], quantile(LRprod[[k]])["50%"], sep=":,")
    report_simu[[k]][[5]] <- paste(simu_items[6], quantile(LRprod[[k]])["75%"], sep=":,")
    report_simu[[k]][[6]] <- paste(simu_items[7], LRprod[[k]][ceiling(length(LRprod[[k]]) * 19 / 20)], sep=":,")
    report_simu[[k]][[7]] <- paste(simu_items[8], quantile(LRprod[[k]])["100%"], sep=":,")
  }

  report_simu.specific <- paste(c("\n\n======== Likelihood ratio ========\n",
                                  "H1 true","\n",
                                  paste(report_simu[[1]], collapse="\n"),
                                  "\n\n","H2 true","\n",
                                  paste(report_simu[[2]], collapse="\n"),
                                  "\n\n======== Parameters ========\n",
                                  "Number of simulations :,",n.simu,"\n",
                                  "\n======== Information ========\n",
                                  "Version: ,", softVer, "\n",
                                  "Date: ,",report_today.simu,"\n",
                                  "Computation time: ,",report_calc.time.simu
  ),sep="",collapse=""
  )
  report <- paste(Report.common("2", dropState, DropMode, FreqMode, FixedLR, PrD, min.freq,
                                freqFp, mutationFp, gtFp,
                                report_locus_set_name, locus_set, df,
                                linkage.loci, RR), report_simu.specific, sep="", collapse="")

  pre_SaveFileName <- tkgetSaveFile(filetypes = "{{CSV files} {.csv}}")
  Report.save(report, pre_SaveFileName)
}


alldata_Report <- function(envProj){
  locus_set <- get("locus_set", pos = envProj)
  linkage.loci <- get("linkage.loci", pos = envProj)
  RR <- get("RR", pos = envProj)
  return_information.simu <- list()
  return_information.simu <- get("return_information.simu", pos = envProj)

  LR.simu <- return_information.simu[[1]]
  calc.time.simu <- return_information.simu[[2]]
  today.simu <- return_information.simu[[3]]

  alldata <- list()
  if (!is.na(RR[1])){
    linkage.loci.vector <- as.vector(linkage.loci)
    alldata.linkage.loci <- apply(linkage.loci, c(1,2), function(x){return(locus_set[x])})
    alldata.allelenames <- paste(c(locus_set[-linkage.loci.vector], paste(alldata.linkage.loci[,1], " - ", alldata.linkage.loci[,2])), sep="", collapse=",")
    alldata.LRprod <- list()
    for (k in 1:2){
      alldata.LRprod[[k]] <- apply(LR.simu[[k]][[1]], 1, prod) * apply(LR.simu[[k]][[2]], 1, prod)
      alldata[[k]] <- paste(paste("H",k,sep=""), apply(LR.simu[[k]][[1]], 1 ,paste, sep="", collapse=","), apply(LR.simu[[k]][[2]], 1, paste, sep="", collapse=","), alldata.LRprod[[k]], sep=",", collapse="\n")
    }
  }else{
    alldata.allelenames <- paste(locus_set, sep="", collapse=",")
    alldata.LRprod <- list()
    for (k in 1:2){
      alldata.LRprod[[k]] <- apply(LR.simu[[k]], 1, prod)
      alldata[[k]] <- paste(paste("H",k,sep=""), apply(LR.simu[[k]], 1, paste, sep="", collapse=","), alldata.LRprod[[k]], sep=",", collapse="\n")
    }
  }

  alldata <- paste(c("True Hypothesis", ",", alldata.allelenames, ",", "Overall LR",
                     "\n",
                     alldata[[1]],
                     "\n",
                     alldata[[2]]
  ), sep="", collapse=""
  )

  pre_SaveFileName <- tkgetSaveFile(filetypes = "{{CSV files} {.csv}}")
  Report.save(alldata, pre_SaveFileName)
}