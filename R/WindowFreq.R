windowFreqSet <- function(envProj){
  FreqMode <- get("FreqMode", pos = envProj)
  FreqModeVar <- tclVar(FreqMode)
  min.freq <- get("min.freq", pos = envProj)
  min.freqVar <- tclVar(min.freq)

  tf <- tktoplevel()
  tkwm.title(tf, "Allele frequency setting")

  outframe_Freq.setting <- tkframe(tf)
  frame_Freq.ok <- tkframe(tf)

  frame_Freq.radiobutton <- tkframe(outframe_Freq.setting)
  frame_Freq.min <- tkframe(outframe_Freq.setting)

  label_min.freq <- tklabel(frame_Freq.min, text="Minimum frequency",state="disable", width="17")
  entry_min.freq <- tkentry(frame_Freq.min, textvariable=min.freqVar, width="15",
                            background="#ffffff" ,state="disable",cursor="arrow")

  if(FreqMode == "2"){
    label_min.freq <- tklabel(frame_Freq.min, text="Minimum frequency",state="normal", width="17")
    entry_min.freq <- tkentry(frame_Freq.min, textvariable=min.freqVar, width="15",
                              background="#ffffff" ,state="normal",cursor="arrow")
  }

  radiobutton_Diri <- tkradiobutton(frame_Freq.radiobutton, anchor="w", width=55, cursor="hand2",
                                    text = "Estimated by Dirichlet distributions", variable = FreqModeVar, value = 0,
                                    command = function() {
                                      tkconfigure(label_min.freq, state="disable")
                                      tkconfigure(entry_min.freq, state="disable")
                                    })
  radiobutton_NRC <- tkradiobutton(frame_Freq.radiobutton, anchor="w", width=55, cursor="hand2",
                                   text = "Corrected by setting the observed number less than five to five", variable = FreqModeVar, value = 1,
                                     command = function() {
                                       tkconfigure(label_min.freq, state="disable")
                                       tkconfigure(entry_min.freq, state="disable")
                                   })
  radiobutton_Act <- tkradiobutton(frame_Freq.radiobutton, anchor="w", width=55, cursor="hand2",
                                   text = "Observed", variable = FreqModeVar, value = 2,
                                   command = function() {
                                               tkconfigure(label_min.freq, state="normal")
                                               tkconfigure(entry_min.freq, state="normal")
                                   })

  tkpack(radiobutton_Diri, anchor="w", padx=15)
  tkpack(radiobutton_NRC, anchor="w", padx=15)
  tkpack(radiobutton_Act, anchor="w", padx=15)
  tkpack(frame_Freq.radiobutton, anchor="w", padx=20, pady=15)

  tkpack(label_min.freq, side="left",anchor="n", padx=10, pady=10)
  tkpack(entry_min.freq, side="left",anchor="n", padx=10, pady=10)
  tkpack(frame_Freq.min, anchor="w", padx=20)

  button_Freq.ok <- tkbutton(frame_Freq.ok,text="Save", width="10", cursor="hand2",
                             command=function() {
                               assign("FreqMode", tclvalue(FreqModeVar), envir = envProj)
                               assign("min.freq", tclvalue(min.freqVar), envir = envProj)
                               tkdestroy(tf)
                             }
  )
  tkpack(button_Freq.ok, padx=15, pady=10)
  tkpack(outframe_Freq.setting, padx=15, pady=10)
  tkpack(frame_Freq.ok, anchor="e", padx=15, pady=10)
}
