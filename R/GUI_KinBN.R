#' KinBN
#'
#' @description Top screen of KinBN
#' @export
KinBN <- function(){

  envProj <- new.env(parent = globalenv())
  setEnvProj(envProj, TRUE)

  envGUI <- new.env(parent = globalenv())

  softVer <- packageVersion("KinBN")
  pathPack <- path.package("KinBN", quiet = FALSE)
  assign("softVer", softVer, envir = envGUI)
  assign("pathPack", pathPack, envir = envGUI)

  linkState <- get("linkState", pos = envProj)
  dropState <- get("dropState", pos = envProj)
  freqState <- get("freqState", pos = envProj)

  linkStateVar <- tclVar(linkState)
  dropStateVar <- tclVar(dropState)
  freqStateVar <- tclVar(freqState)

  modeState <- get("modeState", pos = envProj)
  modeStateVar <- tclVar(modeState)

  tf <- tktoplevel()
  tkwm.title(tf, paste0("KinBN ver. ", softVer))

  topMenu <- tkmenu(tf)
  tkconfigure(tf, menu = topMenu)

  file_menu <- tkmenu(topMenu, tearoff = FALSE, activebackground = "lightskyblue1")
  tkadd(topMenu, "cascade", label = "File", menu = file_menu)
  tkadd(file_menu, "command", label = "New project", command = function() newProj(envProj, envGUI))
  tkadd(file_menu, "command", label = "Load project", command = function() loadProj(tf, envProj, envGUI))
  tkadd(file_menu, "command", label = "Save project", command = function() saveProj(envProj,tclvalue(modeStateVar)))
  tkadd(file_menu, "separator")
  tkadd(file_menu, "command", label = "Quit", command = function() tkdestroy(tf))

  setting_menu <- tkmenu(topMenu, tearoff = FALSE, activebackground = "lightskyblue1")
  tkadd(topMenu, "cascade", label = "Setting", menu = setting_menu)
  tkadd(setting_menu, "command", state = tclvalue(linkStateVar), label = "Linkage", command = function() windowLinkSet(envProj))
  tkadd(setting_menu, "command", state = tclvalue(dropStateVar), label = "Drop-out", command = function() windowDropSet(envProj))
  tkadd(setting_menu, "command", state = tclvalue(freqStateVar), label = "Allele frequency", command = function() windowFreqSet(envProj))
  assign("topMenu", topMenu, envir = envGUI)
  assign("setting_menu", setting_menu, envir = envGUI)

  mode_menu <- tkmenu(topMenu, tearoff = FALSE, activebackground = "lightskyblue1")
  tkadd(topMenu, "cascade", label = "Mode", menu = mode_menu)
  tkadd(mode_menu, "radiobutton", label = "Case analysis", variable = modeStateVar, value = 1, command = function() makeFiles(envProj, envGUI, tclvalue(modeStateVar)))
  tkadd(mode_menu, "radiobutton", label = "Simulation", variable = modeStateVar, value = 2, command = function() makeFiles(envProj, envGUI, tclvalue(modeStateVar)))
  assign("mode_menu", mode_menu, envir = envGUI)

  help_menu <- tkmenu(topMenu, tearoff = FALSE, activebackground = "lightskyblue1")
  tkadd(topMenu, "cascade", label = "Help", menu = help_menu)
  tkadd(help_menu, "command", label = "User manual", command = function() browseURL("https://github.com/ChieMorimoto/KinBN/blob/master/inst/KinBN v2.0.0 user manual.pdf"))
  tkadd(help_menu, "command", label = "About KinBN", command = function() tkmessageBox(title = paste0("KinBN ver. ", softVer), message = "KinBN is a free software GNU General Public License v3.0 for kinship analysis based on a Bayesian network.\nWeb: https://github.com/ChieMorimoto/KinBN", parent = tf, icon="info"))

  tabs <- tk2notebook(tf, tabs = c("Files", "Hypothesis", "Results"))
  tkpack(tabs, fill = "both", expand = 1)
  tabFiles <- tk2notetab(tabs, "Files")
  frameFiles <- tkframe(tabFiles)
  tabHypo <- tk2notetab(tabs, "Hypothesis")
  tabResults <- tk2notetab(tabs, "Results")
  assign("tabs", tabs, envir = envGUI)
  assign("tabFiles", tabFiles, envir = envGUI)
  assign("frameFiles", frameFiles, envir = envGUI)
  assign("tabHypo", tabHypo, envir = envGUI)
  assign("tabResults", tabResults, envir = envGUI)

  frame_Hypo.apply <- tkframe(tabHypo)
  label_ref.num <- tklabel(frame_Hypo.apply, text="Number of known profiles", anchor="w")
  text_ref.num <- tkentry(frame_Hypo.apply, width=8, background="#ffffff")
  button_ref.num <- tkbutton(frame_Hypo.apply, text="Apply", cursor="hand2", width=8)
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


  assign("frame_Hypo.apply", frame_Hypo.apply, envir = envGUI)
  assign("label_ref.num", label_ref.num, envir = envGUI)
  assign("text_ref.num", text_ref.num, envir = envGUI)
  assign("button_ref.num", button_ref.num, envir = envGUI)
  assign("labelframe_Hypo", labelframe_Hypo, envir = envGUI)
  assign("frame_Hypo.button", frame_Hypo.button, envir = envGUI)
  assign("frame_Hypo.label", frame_Hypo.label, envir = envGUI)
  assign("frame_Hypo.next", frame_Hypo.next, envir = envGUI)
  assign("frame_Hypo.select.out",  frame_Hypo.select.out, envir = envGUI)
  assign("canvas_Hypo", canvas_Hypo, envir = envGUI)
  assign("scr_Hypo", scr_Hypo, envir = envGUI)
  assign("frame_Hypo.select", frame_Hypo.select, envir = envGUI)


  labelframe_LR <- tk2labelframe(tabResults,text="Likelihood ratio",labelanchor="nw",relief="groove",borderwidth="2")
  simuframe_LR <- tk2labelframe(tabResults,text="Simulation",labelanchor="nw",relief="groove",borderwidth="2")

  assign("labelframe_LR", labelframe_LR, envir = envGUI)
  assign("simuframe_LR", simuframe_LR, envir = envGUI)

  outframe_simu <- tkframe(tabResults)
  assign("outframe_simu", outframe_simu, envir = envGUI)

  makeFiles(envProj, envGUI, tclvalue(modeStateVar))

}



