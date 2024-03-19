setEnvProj <- function(envProj, setfiles){
  if(setfiles){
    assign("freqInput", NULL, envir = envProj)
    assign("freqFp", character(0), envir = envProj)
    assign("freqFn", character(0), envir = envProj)

    assign("mutationInput", NULL, envir = envProj)
    assign("mutationFp", character(0), envir = envProj)
    assign("mutationFn", character(0), envir = envProj)

    assign("gtInput", NULL, envir = envProj)
    assign("gtFp", character(0), envir = envProj)
    assign("gtFn", character(0), envir = envProj)
  }

  assign("softVer", "2.0.0", envir = envProj)
  assign("ped.max", 25, envir = envProj)

  assign("Mode", "LR", envir = envProj)
  assign("modeState", 1, envir = envProj)

  assign("ref.num", character(0), envir = envProj)

  assign("linkState", "disabled", envir = envProj)
  assign("dropState", "disabled", envir = envProj)
  assign("freqState", "disabled", envir = envProj)
  assign("LinkMode", 1, envir = envProj)
  assign("linkage.loci", NA, envir = envProj)
  assign("RR", NA, envir = envProj)
  assign("linkage_number", 0, envir = envProj)
  assign("locus_set_name", "-", envir = envProj)

  assign("DropMode", 4, envir = envProj)
  assign("FixedLR", 1, envir = envProj)
  assign("PrD", 0.1, envir = envProj)

  assign("simfin", "disabled", envir = envProj)
  assign("FreqMode", 0, envir = envProj)
  assign("min.freq", 0.001, envir = envProj)
  assign("default_locus_set", NULL, envir = envProj)

  assign("freq", NULL, envir = envProj)
  assign("mutation", NULL, envir = envProj)
  assign("gt", NULL, envir = envProj)

  assign("finHypo", FALSE, envir = envProj)
  assign("finResults", FALSE, envir = envProj)
}


newProj <- function(envProj, envGUI){
  inputOk <- "ok"
  finHypo <- get("finHypo", pos = envProj)
  finResults <- get("finResults", pos = envProj)

  if(finHypo || finResults){
    inputOk <- tclvalue(tkmessageBox(message = "Unsaved data will be deleted. Do you want to continue?", type = "okcancel", icon = "warning"))
  }

  if(inputOk == "ok"){
    tabs <- get("tabs", pos = envGUI)

    frameFiles <- get("frameFiles", pos = envGUI)
    frame_Hypo.apply <- get("frame_Hypo.apply", pos = envGUI)
    labelframe_Hypo <- get("labelframe_Hypo", pos = envGUI)
    frame_Hypo.next <- get("frame_Hypo.next", pos = envGUI)
    labelframe_LR <- get("labelframe_LR", pos = envGUI)
    simuframe_LR <- get("simuframe_LR", pos = envGUI)
    outframe_simu <- get("outframe_simu", pos = envGUI)

    tkdestroy(frame_Hypo.apply)
    tkdestroy(labelframe_Hypo[[1]])
    tkdestroy(labelframe_Hypo[[2]])
    tkdestroy(frame_Hypo.next)
    tkdestroy(labelframe_LR)
    tkdestroy(simuframe_LR)
    tkdestroy(outframe_simu)

    envProj <- new.env(parent = globalenv())
    setEnvProj(envProj, TRUE)

    #envGUI <- new.env(parent = globalenv())

    setting_menu <- get("setting_menu", pos = envGUI)
    tkdelete(setting_menu, 0, 2)
    tkadd(setting_menu, "command", state = "disabled", label = "Linkage", command = function() windowLinkSet(envProj))
    tkadd(setting_menu, "command", state = "disabled", label = "Drop-out", command = function() windowDropSet(envProj))
    tkadd(setting_menu, "command", state = "disabled", label = "Allele frequency", command = function() windowFreqSet(envProj))
    assign("linkState", "disabled", envir = envProj)
    assign("dropState", "disabled", envir = envProj)

    modeState <- get("modeState", pos = envProj)
    modeStateVar <- tclVar(modeState)

    mode_menu <- get("mode_menu", pos = envGUI)
    tkdelete(mode_menu, 0, 2)
    tkadd(mode_menu, "radiobutton", label = "Case analysis", variable = modeStateVar, value = 1, command = function() makeFiles(envProj, envGUI, tclvalue(modeStateVar)))
    tkadd(mode_menu, "radiobutton", label = "Simulation", variable = modeStateVar, value = 2, command = function() makeFiles(envProj, envGUI, tclvalue(modeStateVar)))
    makeFiles(envProj, envGUI, tclvalue(modeStateVar))
    tk2notetab.select(tabs, "Files")
  }
}


loadProj <- function(tf, envProj_old, envGUI){
  inputOk <- "ok"
  finHypo <- get("finHypo", pos = envProj_old)
  finResults <- get("finResults", pos = envProj_old)

  if(finHypo || finResults){
    inputOk <- tclvalue(tkmessageBox(message = "Unsaved data will be deleted. Do you want to continue?", type = "okcancel", icon = "warning"))
  }
  if(inputOk == "ok"){
    fpVar <- tclVar("")
    fileName <- tclvalue(tkgetOpenFile(parent = tf, initialdir = tclvalue(fpVar), multiple = "true", filetypes = "{{R Data Files} {.Rdata}}"))
    posPar1 <- regexpr("[{]", fileName)[[1]]
    posPar2 <- regexpr("[}]", fileName)[[1]]
    if((posPar1 == 1) && (posPar2 == nchar(fileName))){
      fileName <- gsub("[{]", "", fileName)
      fileName <- gsub("[}]", "", fileName)
    }
    if(!nchar(fileName)){
      tkmessageBox(message = "No file was selected!", icon = "error", type = "ok")
    }else{
      load(fileName)
      for(i in ls(envProj, all.names = TRUE)){
        assign(i, get(i, envProj), envProj_old)
      }
      modeState_old <- get("modeState", pos = envProj_old)
      modeStateVar_old <- tclVar(modeState_old)
      makeoldFiles(envProj_old, envGUI, tclvalue(modeStateVar_old))
      makeoldHypo(envProj_old, envGUI, tclvalue(modeStateVar_old))
      if(tclvalue(modeStateVar_old) == "1"){
        makeLR(envProj_old, envGUI)
      }else if(tclvalue(modeStateVar_old) == "2"){
        makeSimu(envProj_old, envGUI, tclvalue(modeStateVar_old))
      }
    }
  }
}

# Save a project
saveProj <- function(envProj, modeStateVar){
  assign("modeState", modeStateVar, envir = envProj)
  saveAs <- tkgetSaveFile(filetypes = "{{R Data Files} {.RData}}")
  if(tclvalue(saveAs) != ""){
    if(substr(tclvalue(saveAs), nchar(tclvalue(saveAs)) - 5, nchar(tclvalue(saveAs))) == ".RData"){
      reportName <- tclvalue(saveAs)
    }else{
      reportName <- paste0(tclvalue(saveAs), ".RData")
    }
    save(envProj, file = reportName)
  }
}
