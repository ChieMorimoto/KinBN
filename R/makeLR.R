makeLR <- function(envProj, envGUI){
  tabs <- get("tabs", pos = envGUI)
  tabResults <- get("tabResults", pos = envGUI)
  simfin <- get("simfin", pos = envProj)
  simfinVar <- tclVar(simfin)

  modeState <- get("modeState", pos=envProj)
  modeStateVar <- tclVar(modeState)

  df <- list()
  df[[1]] <- get("df_Hypo1", pos = envProj)
  df[[2]] <- get("df_Hypo2", pos = envProj)
  DNAtype_matrix <- get("DNAtype_matrix", pos = envProj)
  frequency_matrix <- get("frequency_matrix", pos = envProj)
  mutation_matrix <- get("mutation_matrix", pos = envProj)
  locus_set_name <- get("locus_set_name", pos = envProj)
  locus_set <- get("locus_set", pos = envProj)
  linkage.loci <- get("linkage.loci", pos = envProj)
  RR <- get("RR", pos = envProj)
  DropMode <- get("DropMode", pos = envProj)
  FixedLR <- get("FixedLR", pos = envProj)
  PrD <- get("PrD", pos = envProj)
  return_information <- get("return_information", pos = envProj)
  LR <- as.numeric(return_information[[1]])
  calc.time <- return_information[[2]]

  labelframe_LR <- tk2labelframe(tabResults,text="Likelihood ratio",labelanchor="nw",relief="groove",borderwidth="2")
  simuframe_LR <- tk2labelframe(tabResults,text="Simulation",labelanchor="nw",relief="groove",borderwidth="2")

  if (!is.na(RR[1])){
    LRprod <- signif( prod(as.vector(LR)), digits=3)
    log_LR.overall <- paste("Overall LR : ", LRprod)

    linkage.loci.vector <- as.vector(linkage.loci)
    LRtable <- matrix(0, length(locus_set[-linkage.loci.vector]), 2)
    LRtable[,1] <- locus_set[-linkage.loci.vector]
    LRtable[,2] <- signif(LR[-linkage.loci.vector],digits=3)

    LRtable.linkage <- matrix(0, length(RR), 2)
    log_linkage.loci <- apply(linkage.loci,c(1,2),function(x){return(locus_set[x])})
    LRtable.linkage[,1] <- apply(log_linkage.loci, 1 , function(x){return(paste(x[1], " - ", x[2], sep=""))})
    for(i in 1:length(RR)){
      LRtable.linkage[i,2] <- signif(prod(as.numeric(LR[linkage.loci[i,1]]),
                                          as.numeric(LR[linkage.loci[i,2]])),digits=3)
    }
  }else{
    LRprod <- signif(prod(as.vector(LR)),digits=3)
    log_LR.overall <- paste("Overall LR : ", LRprod)
    LRtable <- matrix(0, length(locus_set), 2)
    LRtable[,1] <- locus_set
    LRtable[,2] <- signif(as.vector(LR),digits=3)
    LRtable.linkage <- matrix(0, 1, 1)
    LRtable.linkage[1,1] <- "No linkage"
  }

  label_LR.value <- tklabel(labelframe_LR, text=log_LR.overall)

  box_LRtable <- tk2mclistbox(labelframe_LR, width = 40, height = nrow(LRtable))
  tk2column(box_LRtable, "add", label = "Marker", width = 20)
  tk2column(box_LRtable, "add", label = "LR", width = 20)
  tk2insert.multi(box_LRtable, "end", LRtable)

  box_LRtable.linkage <- tk2mclistbox(labelframe_LR, width = 40, height = nrow(LRtable.linkage))
  tk2column(box_LRtable.linkage, "add", label = "Marker", width = 20)
  tk2column(box_LRtable.linkage, "add", label = "LR", width = 20)
  tk2insert.multi(box_LRtable.linkage, "end", LRtable.linkage)
  button_LR.report <- tkbutton(labelframe_LR, text="Report", width="8", cursor="hand2",
                               command=function() LR_Report(envProj))

  tkpack(label_LR.value, anchor="w", padx=10, pady=10)
  tkpack(box_LRtable, side="left",padx=10, pady=10)
  tkpack(box_LRtable.linkage, side="top",padx=10, pady=10)
  tkpack(button_LR.report, side="bottom", padx=10, pady=10)

  tkpack(labelframe_LR, anchor="w", padx=10, pady=10)

  LR.N <- tclVar("")
  label_LR.N <- tklabel(simuframe_LR, text="Number of simulations",  width="18")
  text_LR.N <- tkentry(simuframe_LR, textvariable=LR.N, width="8",  background="#ffffff")
  button_LR.simu <- tkbutton(simuframe_LR, text="Calculate", cursor="hand2", width="10",
                             command = function() calcSimu(envProj, envGUI, tclvalue(LR.N), tclvalue(modeStateVar)))
  button_LR.simuResult <- tkbutton(simuframe_LR, text="Result", cursor="hand2", width="10", state = tclvalue(simfinVar),
                             command = function() SimuResult(envProj, envGUI))
  tkpack(label_LR.N, side="left",anchor="n", padx=10, pady=10)
  tkpack(text_LR.N, side="left",anchor="n", padx=10, pady=10)
  tkpack(button_LR.simu, side="left",anchor="n", padx=10, pady=10)
  tkpack(button_LR.simuResult, side="left",anchor="n", padx=10, pady=10)
  tkpack(simuframe_LR, anchor="w", padx=10, pady=10)
  assign("LRprod", LRprod, envir = envProj)
  assign("labelframe_LR", labelframe_LR, envir = envGUI)
  assign("simuframe_LR", simuframe_LR, envir = envGUI)
  assign("button_LR.simuResult", button_LR.simuResult, envir = envGUI)
  assign("finResults", TRUE, envir = envProj)
  tk2notetab.select(tabs, "Results")
}

