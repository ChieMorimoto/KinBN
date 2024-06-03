makeFiles <- function(envProj, envGUI, modeStateVar){
    finHypo <- get("finHypo", pos = envProj)
    finResults <- get("finResults", pos = envProj)
    tabs <- get("tabs", pos = envGUI)
    tabHypo <- get("tabHypo", pos = envGUI)
    tabResults <- get("tabResults", pos = envGUI)

    frameFiles <- get("frameFiles", pos = envGUI)
    frame_Hypo.apply <- get("frame_Hypo.apply", pos = envGUI)

    labelframe_Hypo <- get("labelframe_Hypo", pos = envGUI)
    frame_Hypo.next <- get("frame_Hypo.next", pos = envGUI)

    labelframe_LR <- get("labelframe_LR", pos = envGUI)
    simuframe_LR <- get("simuframe_LR", pos = envGUI)

    outframe_simu <- get("outframe_simu", pos = envGUI)

    mode_menu <- get("mode_menu", pos = envGUI)

    tkdestroy(frame_Hypo.apply)

    tkdestroy(labelframe_Hypo[[1]])
    tkdestroy(labelframe_Hypo[[2]])
    tkdestroy(frame_Hypo.next)
    tkdestroy(labelframe_LR)
    tkdestroy(simuframe_LR)
    tkdestroy(outframe_simu)

    setting_menu <- get("setting_menu", pos = envGUI)
    tkdelete(setting_menu, 0, 2)
    tkadd(setting_menu, "command", state = "disabled", label = "Linkage", command = function() windowLinkSet(envProj))
    tkadd(setting_menu, "command", state = "disabled", label = "Drop-out", command = function() windowDropSet(envProj))
    tkadd(setting_menu, "command", state = "disabled", label = "Allele frequency", command = function() windowFreqSet(envProj))
    assign("linkState", "disabled", envir = envProj)
    assign("dropState", "disabled", envir = envProj)
    assign("finHypo", FALSE, envir = envProj)
    assign("finResults", FALSE, envir = envProj)

    tkdestroy(frameFiles)

    freqFn <- get("freqFn", pos = envProj)
    freqFnVar <- tclVar(freqFn)
    mutationFn <- get("mutationFn", pos = envProj)
    mutationFnVar <- tclVar(mutationFn)
    gtFn <- get("gtFn", pos = envProj)
    gtFnVar <- tclVar(gtFn)
    tabFiles <- get("tabFiles", pos = envGUI)

    frameFiles <- tkframe(tabFiles)

    outframe.file <- tkframe(frameFiles)

    tk2notetab.select(tabs, "Files")

    if (modeStateVar == "1") {
      frame_file.DNAtype <- tkframe(outframe.file)
    }
    frame_file.frequency <- tkframe(outframe.file)
    frame_file.mutation <- tkframe(outframe.file)

    frame_Files.next <- tkframe(frameFiles)

    if (modeStateVar == "1"){
      label_DNAtype <- tklabel(frame_file.DNAtype, text="Profiles", width=15,anchor="w")
      check_DNAtype <- tklabel(frame_file.DNAtype, textvariable=gtFnVar,  width=30, relief="groove", borderwidth="3", anchor="w")
      button_DNAtype <- tkbutton(frame_file.DNAtype, text="Browse",  width=10, command=function() OpenFile(envProj, envGUI, "gt", check_DNAtype),cursor="hand2" )
      assign("check_DNAtype", check_DNAtype, envir = envGUI)
    }
    label_frequency <- tklabel(frame_file.frequency, text="Allele frequencies", width=15, anchor="w")
    check_frequency <- tklabel(frame_file.frequency, textvariable=freqFnVar,  width=30, relief="groove", borderwidth="3", anchor="w")
    button_frequency <- tkbutton(frame_file.frequency, text="Browse", width=10, command=function() OpenFile(envProj, envGUI, "freq", check_frequency),cursor="hand2" )
    label_mutation <- tklabel(frame_file.mutation, text="Mutation rates", width=15, anchor="w")
    check_mutation <- tklabel(frame_file.mutation, textvariable=mutationFnVar, width=30, relief="groove", borderwidth="3", anchor="w")
    button_mutation <- tkbutton(frame_file.mutation, text="Browse", width=10, command=function() OpenFile(envProj, envGUI, "mutation", check_mutation),cursor="hand2" )

    button_Files.next <- tkbutton(frame_Files.next,text="Next", width="10",cursor="hand2",
                                  command=function() makeHypo(envProj, envGUI, modeStateVar))

    tkpack(outframe.file,anchor="w",padx="40",pady="20")
    if (modeStateVar == "1"){
      tkpack(frame_file.DNAtype,anchor="w",padx="20",pady="5")
    }
    tkpack(frame_file.frequency,anchor="w",padx="20",pady="5")
    tkpack(frame_file.mutation,anchor="w",padx="20",pady="5")
    tkpack(frame_Files.next, anchor="w", padx="40", pady="5", fill="x")

    if (modeStateVar == "1"){
      tkpack(label_DNAtype,anchor="w",side="left")
      tkpack(check_DNAtype,anchor="w",side="left",padx="5")
      tkpack(button_DNAtype,padx="5")
    }
    tkpack(label_frequency,anchor="w",side="left")
    tkpack(check_frequency,anchor="w",side="left",padx="5")
    tkpack(button_frequency,padx="5")
    tkpack(label_mutation,anchor="w",side="left")
    tkpack(check_mutation,anchor="w",side="left",padx="5")
    tkpack(button_mutation,padx="5")

    tkpack(button_Files.next, anchor="e", padx=10, pady=10)

    tkgrid(frameFiles)
    assign("frameFiles", frameFiles, envir = envGUI)

    assign("check_frequency", check_frequency, envir = envGUI)
    assign("check_mutation", check_mutation, envir = envGUI)
}

