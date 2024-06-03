Report.common <- function(modeStateVar, dropState, DropMode, FreqMode, FixedLR, PrD, min.freq, freqFp, mutationFp, gtFp,
                          report_locus_set_name, locus_set, df,
                          linkage.loci, RR){

  report_locus_set <- paste(locus_set, sep="", collapse=",")
  result_pedigree.array <- list()
  for (j in 1:2){
    result_pedigree.array[[j]] <- array("-", dim=c(nrow(df[[j]]), 3))
    for (i in 1:nrow(df[[j]])){
      result_pedigree.array[[j]][i, 1] <- as.character(df[[j]]$names)[i]
      if (df[[j]]$pat[i] != 0){
        result_pedigree.array[[j]][i, 2] <- as.character(df[[j]]$names)[df[[j]]$pat[i]]
      }
      if (df[[j]]$mat[i] != 0){
        result_pedigree.array[[j]][i, 3] <- as.character(df[[j]]$names)[df[[j]]$mat[i]]
      }
    }
  }
  report_pedigree <- paste("H1\nName,Father,Mother\n",
                           paste(apply(result_pedigree.array[[1]], 1, function(x){paste(x[1], x[2], x[3], sep=",")}), sep="", collapse="\n"),"\n",
                           "\nH2\nName,Father,Mother\n",
                           paste(apply(result_pedigree.array[[2]], 1, function(x){paste(x[1], x[2], x[3], sep=",")}), sep="", collapse="\n"), sep=""
  )

  if (!is.na(RR[1])){
    report.linkage.loci <- apply(linkage.loci,c(1,2),function(x){return(locus_set[x])})
    report.linkage.allelenames <- paste("Locus 1,Locus 2,Recombination rate\n",
                                        paste(paste(apply(report.linkage.loci, 1, paste, collapse=","), RR, sep=","), collapse="\n"), sep="", collapse="")
  }else{
    report.linkage.allelenames <- "No linkage"
  }

  if (DropMode == "4"){
    report.drop <- "No drop-out"
  }else if(DropMode == "0"){
    report.drop <- paste("User-defied LR", "(LR=", FixedLR, ")", sep="")
  }else if(DropMode == "2"){
    report.drop <- "LR considering all possible genotypes equally (Method A)"
  }else if(DropMode == "3"){
    report.drop <- paste("LR calculated with Pr(D) (Method B)", "(Pr(D)=", PrD, ")", sep="")
  }
  if (FreqMode=="0"){
    report.allelefreq <- "Estimated by Dirichlet distributions"
  }else if(FreqMode=="1"){
    report.allelefreq <- "Corrected by setting the observed number less than five to five"
  }else if(FreqMode=="2"){
    report.allelefreq <- paste("Observed", "(minimum frequency=", min.freq, ")", sep="")
  }
  if (modeStateVar == "1"){
    report_common <- paste(c("======== Files ========\n",
                             "Allele frequencies,",freqFp,"\n",
                             "Mutation rates,",mutationFp,"\n",
                             "Profiles,",gtFp,"\n",
                             "\n======== Locus set ========\n",
                             "Kit,",report_locus_set_name,"\n",
                             "Locus,",report_locus_set,"\n",
                             "\n======== Linkage ========\n",
                             report.linkage.allelenames,"\n",
                             "\n======== Drop-out ========\n",
                             report.drop,"\n",
                             "\n======== Allele frequency ========\n",
                             report.allelefreq,"\n",
                             "\n======== Hypothesis ========\n",
                             report_pedigree
    ),sep="",collapse=""
    )
  }else{
    report_common <- paste(c("======== Files ========\n",
                             "Allele frequencies,",freqFp,"\n",
                             "Mutation rates,",mutationFp,"\n",
                             "\n======== Locus set ========\n",
                             "Kit,",report_locus_set_name,"\n",
                             "Locus,",report_locus_set,"\n",
                             "\n======== Linkage ========\n",
                             report.linkage.allelenames,"\n",
                             "\n======== Allele frequency ========\n",
                             report.allelefreq,"\n",
                             "\n======== Hypothesis ========\n",
                             report_pedigree
    ),sep="",collapse=""
    )
  }
  return(report_common)
}

Report.save <- function(report, pre_SaveFileName){
  if (tclvalue(pre_SaveFileName) != ""){
    if (substr(tclvalue(pre_SaveFileName),nchar(tclvalue(pre_SaveFileName)) - 3,nchar(tclvalue(pre_SaveFileName))) == ".csv"){
      SaveFileName <- tclvalue(pre_SaveFileName)
      if (file.exists(SaveFileName)){
        write(report,file=SaveFileName)
      }else{
        file.create(SaveFileName)
        write(report,file=SaveFileName)
      }
    }else{
      SaveFileName <- paste(tclvalue(pre_SaveFileName),".csv",sep="")
      if (file.exists(SaveFileName)){
        warningMessage <- paste(strsplit(SaveFileName,"/")[[1]][length(strsplit(SaveFileName,"/")[[1]])]," already exists.\nDo you want to overwrite it?")
        overwrite <- tkmessageBox(message=warningMessage, type="okcancel", icon="warning")
        if (tclvalue(overwrite) == "ok"){
          write(report,file=SaveFileName)
        }
      }else{
        file.create(SaveFileName)
        write(report,file=SaveFileName)
      }
    }
  }
}

LR_Report <- function(envProj){
  #modeState <- get("modeState", pos = envProj)
  #modeStateVar <- tclVar(modeState)
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
  return_information <- list()
  return_information <- get("return_information", pos = envProj)

  LR <- as.vector(as.numeric(return_information[[1]]))
  calc.time <- return_information[[2]]
  today <- return_information[[3]]
  df <- list()
  df[[1]] <- get("df_Hypo1", pos = envProj)
  df[[2]] <- get("df_Hypo2", pos = envProj)

  if (!is.na(RR[1])){
    linkage.loci.vector <- as.vector(linkage.loci)
    report_LRprod <- prod(LR)
    report_result <- paste(locus_set[-linkage.loci.vector], ",", LR[-linkage.loci.vector], sep="", collapse="\n")

    report_linkage.loci <- apply(linkage.loci,c(1,2),function(x){return(locus_set[x])})
    report_linkage_each <- rep(0,length(RR))
    for(i in 1:length(RR)){
      report_linkage_each[i] <- paste(report_linkage.loci[i,1], " - ", report_linkage.loci[i,2], ",",
                                    prod(LR[linkage.loci[i,1]], LR[linkage.loci[i,2]]))
      report_linkage <- paste(report_linkage_each, collapse="\n")
    }

  }else{
    report_result <- paste(paste(locus_set, ",", LR, sep=""), sep="", collapse="\n")
    report_LRprod <- prod(LR)
    report_linkage <- ""
  }
  report_today <- format(today, "%Y/%b/%d %X")
  report_calc.time <- paste(signif(calc.time[3], digits=3), "sec", collapse="")

  report_LR.specific <- paste(c("\n\n======== Likelihood ratio ========\n",
                                "Overall LR, ",report_LRprod,"\n\n",
                                report_result,"\n\n",
                                report_linkage,"\n",
                                "\n======== Information ========\n",
                                "Version,", softVer, "\n",
                                "Date,",report_today,"\n",
                                "Computation time,",report_calc.time
  ),sep="",collapse=""
  )

  report <- paste(Report.common("1", dropState, DropMode, FreqMode, FixedLR, PrD, min.freq, freqFp, mutationFp, gtFp,
                                report_locus_set_name, locus_set, df,
                                linkage.loci, RR), report_LR.specific, sep="", collapse="")

  pre_SaveFileName <- tkgetSaveFile(filetypes = "{{CSV files} {.csv}}")
  Report.save(report, pre_SaveFileName)
}
