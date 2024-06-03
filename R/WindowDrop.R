windowDropSet <- function(envProj){
  locus_set <- get("locus_set", pos = envProj)
  locus_set_name <- get("locus_set_name", pos = envProj)
  pair.max <- floor(length(locus_set) / 2)
  DropMode <- get("DropMode", pos = envProj)
  DropModeVar <- tclVar(DropMode)

  assign("DropMode", DropMode, envir = envProj)
  FixedLR <- get("FixedLR", pos = envProj)
  FixedLRVar <- tclVar(FixedLR)
  PrD <- get("PrD", pos = envProj)
  PrDVar <- tclVar(PrD)

  tf <- tktoplevel()
  tkwm.title(tf, "Drop-out setting")

  outframe_Drop.setting <- tkframe(tf)
  frame_Drop.ok <- tkframe(tf)

  frame_Drop.radiobutton <- tkframe(outframe_Drop.setting)

  frame_FixedLR <- tkframe(outframe_Drop.setting)
  frame_Semi <- tkframe(outframe_Drop.setting)

  label_FixedLR <- tklabel(frame_FixedLR, text="LR value in allele drop-out locus", width="30",state="disabled")
  entry_FixedLR <- tkentry(frame_FixedLR, textvariable = FixedLRVar, width="5", bg="white",state="disabled")
  if(DropMode == "0"){
    label_FixedLR <- tklabel(frame_FixedLR, text="LR value in allele drop-out locus", width="30",state="normal")
    entry_FixedLR <- tkentry(frame_FixedLR, textvariable=FixedLRVar, width="5",
                              background="#ffffff" ,state="normal",cursor="arrow")
  }

  label_Semi <- tklabel(frame_Semi,text="Pr(D)", width="8",state="disabled")
  entry_Semi <- tkentry(frame_Semi, textvariable = PrDVar, width="5", bg="white",state="disabled")
  if(DropMode == "3"){
    label_Semi <- tklabel(frame_Semi,text="Pr(D)", width="8",state="normal")
    entry_Semi <- tkentry(frame_Semi, textvariable = PrDVar, width="5",
                             background="#ffffff" ,state="normal",cursor="arrow")
  }

  radiobutton_FixLR <- tkradiobutton(frame_Drop.radiobutton, anchor="w", width=50, cursor="hand2",
                                     text = "User-defined LR", variable = DropModeVar, value = 0,
                                     command = function() {
                                       tkconfigure(label_FixedLR, state="normal")
                                       tkconfigure(entry_FixedLR, state="normal")
                                       tkconfigure(label_Semi, state="disabled")
                                       tkconfigure(entry_Semi, state="disabled")
                                     })
  radiobutton_Both <- tkradiobutton(frame_Drop.radiobutton, anchor="w", width=50, cursor="hand2",
                                    text = "LR considering all possible genotypes equally (Method A)", variable = DropModeVar, value = 2,
                                    command = function() {
                                      tkconfigure(label_FixedLR, state="disabled")
                                      tkconfigure(entry_FixedLR, state="disabled")
                                      tkconfigure(label_Semi, state="disabled")
                                      tkconfigure(entry_Semi, state="disabled")
                                    })
  radiobutton_Semi <- tkradiobutton(frame_Drop.radiobutton, anchor="w", width=50, cursor="hand2",
                                    text = "LR calculated with Pr(D) (Method B)", variable = DropModeVar, value = 3,
                                    command = function() {
                                      tkconfigure(label_FixedLR, state="disabled")
                                      tkconfigure(entry_FixedLR, state="disabled")
                                      tkconfigure(label_Semi, state="normal")
                                      tkconfigure(entry_Semi, state="normal")
                                    })
  tkpack(radiobutton_Both, anchor="w", padx=15)

  tkpack(radiobutton_Semi, anchor="w", padx=15)
  tkpack(radiobutton_FixLR, anchor="w", padx=15)

  tkpack(frame_Drop.radiobutton, anchor="e", padx=20, pady=15)

  tkpack(label_Semi, side="left",anchor="w", padx=10, pady=10)
  tkpack(entry_Semi, side="left",anchor="n", padx=10, pady=10)
  tkpack(frame_Semi, anchor="w", padx=15)

  tkpack(label_FixedLR, side="left", anchor="w", padx=10, pady=10)
  tkpack(entry_FixedLR, side="left",anchor="n", padx=10, pady=10)
  tkpack(frame_FixedLR, anchor="w", padx=15)

  button_Drop.ok <- tkbutton(frame_Drop.ok,text="Save", width="8", cursor="hand2",
                             command=function() {
                               assign("DropMode", tclvalue(DropModeVar), envir = envProj)
                               assign("FixedLR", tclvalue(FixedLRVar), envir = envProj)
                               assign("PrD", tclvalue(PrDVar), envir = envProj)
                               tkdestroy(tf)
                             }
  )
  tkpack(button_Drop.ok, anchor="e", padx=15, pady=10)
  tkpack(outframe_Drop.setting, padx=15, pady=10)
  tkpack(frame_Drop.ok, anchor="e", padx=15, pady=10)
}



