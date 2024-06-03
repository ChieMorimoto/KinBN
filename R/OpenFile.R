OpenFile <- function(envProj, envGUI, filetype, filelabel){
  inputOk <- "ok"
  finHypo <- get("finHypo", pos = envProj)
  finResults <- get("finResults", pos = envProj)
  labelframe_Hypo <- get("labelframe_Hypo", pos = envGUI)
  frame_Hypo.button <- get("frame_Hypo.button", pos = envGUI)
  frame_Hypo.label <- get("frame_Hypo.label", pos = envGUI)
  frame_Hypo.next <- get("frame_Hypo.next", pos = envGUI)
  labelframe_LR <- get("labelframe_LR", pos = envGUI)
  simuframe_LR <- get("simuframe_LR", pos = envGUI)
  linkState <- get("linkState", pos = envProj)
  dropState <- get("dropState", pos = envProj)
  linkStateVar <- tclVar(linkState)
  dropStateVar <- tclVar(dropState)
  modeState <- get("modeState", pos = envProj)
  modeStateVar <- tclVar(modeState)
  setting_menu <- get("setting_menu", pos = envGUI)

  if(finHypo){
    if(finHypo || finResults){
      inputOk <- tclvalue(tkmessageBox(message = "Setting of Hypothesis will be deleted. Do you want to continue?", type = "okcancel", icon = "warning"))
    }
    if(inputOk == "ok"){
      tclvalue(linkStateVar) <-  "disabled"
      assign("linkState", tclvalue(linkStateVar), envir = envProj)
      tclvalue(dropStateVar) <-  "disabled"
      assign("dropState", tclvalue(dropStateVar), envir = envProj)
      tkdelete(setting_menu, 0, 2)
      tkadd(setting_menu, "command", state = tclvalue(linkStateVar), label = "Linkage", command = function() windowLinkSet(envProj))
      tkadd(setting_menu, "command", state = tclvalue(dropStateVar), label = "Drop-out", command = function() windowDropSet(envProj))
      tkadd(setting_menu, "command", state = "disabled", label = "Allele frequency", command = function() windowFreqSet(envProj))
      for (j in 1:2){
        tkdestroy(labelframe_Hypo[[j]])
        tkdestroy(frame_Hypo.button[[j]])
        tkdestroy(frame_Hypo.label[[j]])
      }
      tkdestroy(frame_Hypo.next)
      tkdestroy(labelframe_LR)
      tkdestroy(simuframe_LR)

      tabHypo <- get("tabHypo", pos = envGUI)
      tabResults <- get("tabResults", pos = envGUI)
      labelframe_Hypo <- list()
      frame_Hypo.button <- list()
      frame_Hypo.label <- list()
      for (j in 1:2){
        labelframe_Hypo[[j]] <- tk2labelframe(tabHypo,text=paste("Hypothesis", j, collapse="", sep=""),labelanchor="nw",relief="groove",borderwidth="2")
        frame_Hypo.button[[j]] <- tkframe(labelframe_Hypo[[j]])
        frame_Hypo.label[[j]] <- tkframe(labelframe_Hypo[[j]])
      }
      frame_Hypo.next <- tkframe(tabHypo)

      frame_Hypo.select.out <- vector("list", length = 2)
      canvas_Hypo <- vector("list", length = 2)
      scr_Hypo <- vector("list", length = 2)
      frame_Hypo.select <- vector("list", length = 2)
      for (j in 1:2){
        frame_Hypo.select.out[[j]] <- tkframe(labelframe_Hypo[[j]])
        canvas_Hypo[[j]] <- tkcanvas(frame_Hypo.select.out[[j]])
        scr_Hypo[[j]] <- tkscrollbar(frame_Hypo.select.out[[j]])
        frame_Hypo.select[[j]] <- tkframe(canvas_Hypo[[j]])
      }

      assign("labelframe_Hypo", labelframe_Hypo, envir = envGUI)
      assign("frame_Hypo.button", frame_Hypo.button, envir = envGUI)
      assign("frame_Hypo.label", frame_Hypo.label, envir = envGUI)
      assign("frame_Hypo.next", frame_Hypo.next, envir = envGUI)
      assign("frame_Hypo.select.out",  frame_Hypo.select.out, envir = envGUI)
      assign("canvas_Hypo", canvas_Hypo, envir = envGUI)
      assign("scr_Hypo", scr_Hypo, envir = envGUI)
      assign("frame_Hypo.select", frame_Hypo.select, envir = envGUI)


      labelframe_LR <- tk2labelframe(tabResults,text="Likelihood ratio",labelanchor="nw",relief="groove",borderwidth="2")
      frame_LR.locus.value <- tkframe(labelframe_LR)
      frame_LR.linkage <- tkframe(labelframe_LR)
      frame_LR.linkage.locus.value <- tkframe(frame_LR.linkage)
      frame_LR.report <- tkframe(frame_LR.linkage)
      assign("labelframe_LR", labelframe_LR, envir = envGUI)
      assign("frame_LR.locus.value", frame_LR.locus.value, envir = envGUI)
      assign("frame_LR.linkage", frame_LR.linkage, envir = envGUI)
      assign("frame_LR.linkage.locus.value", frame_LR.linkage.locus.value, envir = envGUI)
      assign("frame_LR.report", frame_LR.report, envir = envGUI)

      setEnvProj(envProj, FALSE)
      makeFiles(envProj, envGUI, tclvalue(modeStateVar))
    }
  }

  if(inputOk == "ok"){
    if(filetype == "freq"){
      fp <- get("freqFp", pos = envProj)
      fn <- get("freqFn", pos = envProj)
    }else if(filetype == "mutation"){
      fp <- get("mutationFp", pos = envProj)
      fn <- get("mutationFn", pos = envProj)
    }else if(filetype == "gt"){
      fp <- get("gtFp", pos = envProj)
      fn <- get("gtFn", pos = envProj)
    }
    fpVar <- tclVar(fp)
    fnVar <- tclVar(fn)

    fileName <-  tclvalue(tkgetOpenFile(initialdir = tclvalue(fpVar), multiple = "true", filetypes = "{{CSV Files} {.csv}}"))
    if(!nchar(fileName)){
      tkmessageBox(message = "No file was selected!", icon = "error", type = "ok")
    }else{
      tmp <- sub("\\}", fileName, replacement = "")
      tmp2 <- sub("\\{", tmp, replacement = "")
      tclvalue(fpVar) <- tmp2
      foo3 <- strsplit(tmp2, "/")[[1]]
      tclvalue(fnVar) <- strsplit(foo3[length(foo3)], "\\.csv")[[1]][1]

      if(filetype == "freq"){
        freqInput <- read.csv(tclvalue(fpVar), header = TRUE)
        freqInput <- as.matrix(freqInput)
        freqInput <- gsub(" ", "", freqInput, fixed = TRUE)
        assign("freqInput", freqInput, envir = envProj)
        assign("freqFp", tclvalue(fpVar), envir = envProj)
        assign("freqFn", tclvalue(fnVar), envir = envProj)
        filelabel <- get("check_frequency", pos = envGUI)
      }else if(filetype == "mutation"){
        mutationInput <- read.csv(tclvalue(fpVar), header = TRUE)
        mutationInput <- as.matrix(mutationInput)
        mutationInput <- gsub(" ", "", mutationInput, fixed = TRUE)
        assign("mutationInput", mutationInput, envir = envProj)
        assign("mutationFp", tclvalue(fpVar), envir = envProj)
        assign("mutationFn", tclvalue(fnVar), envir = envProj)
        filelabel <- get("check_mutation", pos = envGUI)
      }else if(filetype == "gt"){
        gtInput <- read.csv(tclvalue(fpVar), header = TRUE)
        gtInput <- as.matrix(gtInput)
        colnames(gtInput) <- gsub(".", "", colnames(gtInput), fixed = TRUE)
        assign("gtInput", gtInput, envir = envProj)
        assign("gtFp", tclvalue(fpVar), envir = envProj)
        assign("gtFn", tclvalue(fnVar), envir = envProj)
        filelabel <- get("check_DNAtype", pos = envGUI)
      }
      tkconfigure(filelabel, textvariable = fnVar)
    }
  }

}
