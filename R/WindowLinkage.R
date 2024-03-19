windowLinkSet <- function(envProj){
  LinkMode <- get("LinkMode", pos = envProj)
  LinkModeVar <- tclVar(LinkMode)
  linkage.loci <- get("linkage.loci", pos = envProj)
  RR <- get("RR", pos = envProj)
  linkage_number <- get("linkage_number", pos = envProj)
  var_locus <- locus_set <- get("locus_set", pos = envProj)
  locus_set_name <- get("locus_set_name", pos = envProj)
  pair.max <- floor(length(locus_set) / 2)
  default_locus_set <- get("default_locus_set", pos = envProj)
  locus1.var <- locus2.var <- RR.var <- list()
  combobox_Linkage.locus1 <- combobox_Linkage.locus2 <- entry_Linkage.RR <- button_Linkage.vanish <- list()

  button_vanish <- function(i){
    force(i);function(){
      if (i != linkage_number){
        for (k in (i + 1):linkage_number){
          tclvalue(locus1.var[[k - 1]]) <<- tclvalue(locus1.var[[k]])
          tclvalue(locus2.var[[k - 1]]) <<- tclvalue(locus2.var[[k]])
          tclvalue(RR.var[[k - 1]]) <<- tclvalue(RR.var[[k]])
        }
      }
      tkpack.forget(frame_Linkage.select.input[[linkage_number]])
      tclvalue(locus1.var[[linkage_number]]) <<- ""
      tclvalue(locus2.var[[linkage_number]]) <<- ""
      tclvalue(RR.var[[linkage_number]]) <<- ""
      linkage_number <<- linkage_number - 1
    }
  }

  cb_default.change <- function(x){
    if (x == "1"){
      tkpack(frame_Linkage.add, anchor="w", padx=40, pady=5)
      tkpack(frame_Linkage.label, anchor="w", padx=40, pady=5)
      tkpack(frame_Linkage.select, padx=40, anchor="w")
      if (linkage_number > 0){
        for (i in 1:linkage_number){
          tkpack.forget(frame_Linkage.select.input[[i]])
        }
      }
      linkage_number.default <<- length(RR.default)
      default.sort <- sapply(default_locus_set[[locus_set_name]][[1]], function(x){return(which(x == locus_set))})
      for (i in 1:linkage_number.default){
        tkpack(frame_Linkage.select.input[[i]], anchor="w", pady="5")
        tkconfigure(combobox_Linkage.locus1[[i]], state="disable")
        tkconfigure(combobox_Linkage.locus2[[i]], state="disable")
        tkconfigure(entry_Linkage.RR[[i]], state="disable",cursor="arrow")
        tkconfigure(button_Linkage.vanish[[i]], state="disable",cursor="arrow")
        tclvalue(locus1.var[[i]]) <- var_locus[default.sort[linkage.loci.default[i,1]]]
        tclvalue(locus2.var[[i]]) <- var_locus[default.sort[linkage.loci.default[i,2]]]
        tclvalue(RR.var[[i]]) <- RR.default[i]
      }
      tkconfigure(button_Linkage.add, state="disable",cursor="arrow")
    }else if(x == "2"){
      tkpack(frame_Linkage.add, anchor="w", padx=40, pady=5)
      tkpack(frame_Linkage.label, anchor="w", padx=40, pady=5)
      tkpack(frame_Linkage.select, padx=40, anchor="w")
      for (i in 1:pair.max){
        tkconfigure(combobox_Linkage.locus1[[i]], state="normal")
        tkconfigure(combobox_Linkage.locus2[[i]], state="normal")
        tkconfigure(entry_Linkage.RR[[i]], state="normal",cursor="xterm")
        tkconfigure(button_Linkage.vanish[[i]], state="normal",cursor="hand2")
      }
      tkconfigure(button_Linkage.add, state="normal",cursor="hand2")
    }else{
      tkpack.forget(frame_Linkage.add)
      tkpack.forget(frame_Linkage.label)
      tkpack.forget(frame_Linkage.select)
    }
  }


  checkLinkage <- function(locus_set_name, LinkModeVar, var_locus, linkage_number, locus1.var, locus2.var, RR.var){
    is.default <- sapply(default_locus_set, function(x){if (setequal(x[[1]], locus_set) && length(x[[1]]) == length(locus_set) ) return(1) else return(0)})

    if (LinkModeVar == "1"){
      locus_set_name <- names(default_locus_set)[which(is.default == 1)]
      default.sort <- sapply(default_locus_set[[locus_set_name]][[1]], function(x){return(which(x == locus_set))})
      linkage.loci <- apply(default_locus_set[[locus_set_name]][[2]], c(1,2), function(x){return(default.sort[x])})
      RR <- default_locus_set[[locus_set_name]][[3]]
      return(list(linkage.loci, RR))
    }else if (LinkModeVar == "3"){
      linkage.loci <- RR <- NA
      #locus_set_name <- "-"
      return(list(linkage.loci, RR))
    }else if(LinkModeVar == "2"){
      if(linkage_number == 0){
        linkage_number <- 1
      }
      linkage.pair <- sapply(1:linkage_number,function(x){if(tclvalue(locus1.var[[x]]) != "" &&
                                                             tclvalue(locus2.var[[x]]) != "" &&
                                                             tclvalue(RR.var[[x]]) != "") return(x) else return(0)})
     if (is.element(0, linkage.pair)){
       tkmessageBox(message = "Error: incorrect setting of linked loci", icon = "error")
      }else{
        locus1 <- sapply(linkage.pair,function(x){return(which(var_locus == tclvalue(locus1.var[[x]])))})
        locus2 <- sapply(linkage.pair,function(x){return(which(var_locus == tclvalue(locus2.var[[x]])))})
        RR <- sapply(linkage.pair,function(x){return(as.numeric(tclvalue(RR.var[[x]])))})

        if (!length(unique(c(locus1,locus2))) == length(c(locus1,locus2))){
          tkmessageBox(message = "Error2", icon = "error")

          tkmessageBox(message = "Error: incorrect setting of linked loci", icon="error")
        }else if (prod(0 <= RR) * prod(RR <= 1) == 0){
          tkmessageBox(message = "Please enter a value of 0-1 to Recombination rate.", icon="error")
        }else{
          linkage.loci <- cbind(locus1,locus2)
          #locus_set_name <- "-"
          return(list(linkage.loci, RR))
        }
      }
    }
  }


  tf <- tktoplevel()
  tkwm.title(tf, "Linkage setting")

  outframe_Linkage.setting <- tkframe(tf)

  frame_Linkage.radiobutton <- tkframe(outframe_Linkage.setting)
  frame_Linkage.add <- tkframe(outframe_Linkage.setting)
  frame_Linkage.label <- tkframe(outframe_Linkage.setting)
  frame_Linkage.select <- tkframe(outframe_Linkage.setting)

  frame_Linkage.select.input <- list()

  button_Linkage.add <- tkbutton(frame_Linkage.add,text="Add",width=10,state="disable",cursor="arrow",
                                 command = function(){
                                   if (0 <= linkage_number && linkage_number < pair.max){
                                     linkage_number <<- linkage_number + 1
                                     tkconfigure(combobox_Linkage.locus1[[linkage_number]], state="normal")
                                     tkconfigure(combobox_Linkage.locus2[[linkage_number]], state="normal")
                                     tkconfigure(entry_Linkage.RR[[linkage_number]], state="normal")
                                     tkconfigure(button_Linkage.vanish[[linkage_number]], state="normal")
                                     tkpack(frame_Linkage.select.input[[linkage_number]],anchor="w", pady="5")
                                   }
                                 }
  )

  for (i in 1:pair.max){
    locus1.var[[i]] <- tclVar("")
    locus2.var[[i]] <- tclVar("")
    RR.var[[i]] <- tclVar("")
    frame_Linkage.select.input[[i]] <- tkframe(frame_Linkage.select)
    combobox_Linkage.locus1[[i]] <- ttkcombobox(frame_Linkage.select.input[[i]],values=var_locus,textvariable=locus1.var[[i]],width="10",state="disable")
    combobox_Linkage.locus2[[i]] <- ttkcombobox(frame_Linkage.select.input[[i]],values=var_locus,textvariable=locus2.var[[i]],width="10",state="disable")
    entry_Linkage.RR[[i]] <- tkentry(frame_Linkage.select.input[[i]], textvariable=RR.var[[i]], width="15",  background="#ffffff" ,state="disable",cursor="arrow")
    button_Linkage.vanish[[i]] <- tkbutton(frame_Linkage.select.input[[i]],text="x",width=5,state="disable", command = button_vanish(i))
  }

  if (locus_set_name != "-"){
    radiobutton_off <- tkradiobutton(frame_Linkage.radiobutton, anchor="w", width=40, cursor="hand2",
                                     text=paste("Default settings of", locus_set_name, collapse=""), variable = LinkModeVar,
                                     value = 1, command=function(){cb_default.change(as.integer(tclvalue(LinkModeVar)))})
    radiobutton_on <- tkradiobutton(frame_Linkage.radiobutton, anchor="w", width=40, cursor="hand2",
                                    text="Manual setting", variable = LinkModeVar,
                                    value = 2, command=function(){cb_default.change(as.integer(tclvalue(LinkModeVar)))})
    radiobutton_none <- tkradiobutton(frame_Linkage.radiobutton, anchor="w", width=40, cursor="hand2",
                                      text="No linkage", variable = LinkModeVar,
                                      value = 3, command=function(){cb_default.change(as.integer(tclvalue(LinkModeVar)))})

    linkage.loci.default <- default_locus_set[[locus_set_name]][[2]]
    RR.default <- default_locus_set[[locus_set_name]][[3]]
    linkage_number.default <- length(RR.default)
    default.sort <- sapply(default_locus_set[[locus_set_name]][[1]], function(x){return(which(x == locus_set))})

    if(tclvalue(LinkModeVar) == "2"){
      for (i in 1:nrow(linkage.loci)){
        linkage_number <- nrow(linkage.loci)
        tkpack(frame_Linkage.select.input[[i]], anchor="w", pady=5)
        tclvalue(locus1.var[[i]]) <- var_locus[linkage.loci[i,1]]
        tclvalue(locus2.var[[i]]) <- var_locus[linkage.loci[i,2]]
        tclvalue(RR.var[[i]]) <- RR[i]
        tkconfigure(combobox_Linkage.locus1[[i]], state="normal")
        tkconfigure(combobox_Linkage.locus2[[i]], state="normal")
        tkconfigure(entry_Linkage.RR[[i]], state="normal")
      }
    }else{
      for (i in 1:linkage_number.default){
        tkpack(frame_Linkage.select.input[[i]],anchor="w", pady=5)
        tclvalue(locus1.var[[i]]) <- var_locus[default.sort[linkage.loci.default[i,1]]]
        tclvalue(locus2.var[[i]]) <- var_locus[default.sort[linkage.loci.default[i,2]]]
        tclvalue(RR.var[[i]]) <- RR.default[i]
      }
    }


  }else if(locus_set_name == "-"){
    if(tclvalue(LinkModeVar) == "1"){
      LinkModeVar <- tclVar("2")
    }
    radiobutton_on <- tkradiobutton(frame_Linkage.radiobutton, anchor="w", width=40, cursor="hand2",
                                    text="Manual setting", variable=LinkModeVar, value=2, command=function(){cb_default.change(as.integer(tclvalue(LinkModeVar)))})
    radiobutton_none <- tkradiobutton(frame_Linkage.radiobutton, anchor="w", width=40, cursor="hand2",
                                      text="No linkage", variable=LinkModeVar, value=3, command=function(){cb_default.change(as.integer(tclvalue(LinkModeVar)))})
    if(!is.na(RR[1])){
      for (i in 1:nrow(linkage.loci)){
        linkage_number <- nrow(linkage.loci)
        tkpack(frame_Linkage.select.input[[i]], anchor="w", pady=5)
        tclvalue(locus1.var[[i]]) <- var_locus[linkage.loci[i,1]]
        tclvalue(locus2.var[[i]]) <- var_locus[linkage.loci[i,2]]
        tclvalue(RR.var[[i]]) <- RR[i]
        tkconfigure(combobox_Linkage.locus1[[i]], state="normal")
        tkconfigure(combobox_Linkage.locus2[[i]], state="normal")
        tkconfigure(entry_Linkage.RR[[i]], state="normal")
        tkpack(frame_Linkage.select.input[[i]], anchor="w", pady=5)
      }

    }else{
      for (i in 1:pair.max){
        tkconfigure(combobox_Linkage.locus1[[i]], state="readonly")
        tkconfigure(combobox_Linkage.locus2[[i]], state="readonly")
        tkconfigure(entry_Linkage.RR[[i]], state="normal")
      }
      tkpack(frame_Linkage.select.input[[1]], anchor="w", pady=5)
    }

    tkconfigure(button_Linkage.add, state="normal",cursor="hand2")
  }

  frame_Linkage.save <- tkframe(tf)
  button_Linkage.save <- tkbutton(frame_Linkage.save,text="Save",width="10",cursor="hand2",
                                  command = function(){
                                    LinkageSetting <- checkLinkage(locus_set_name, tclvalue(LinkModeVar),
                                                                    var_locus, linkage_number,
                                                                    locus1.var, locus2.var,
                                                                    RR.var)
                                    linkage.loci <- LinkageSetting[[1]]
                                    RR <- LinkageSetting[[2]]
                                    assign("linkage.loci", linkage.loci, envir = envProj)
                                    assign("RR", RR, envir = envProj)
                                    assign("linkage_number", linkage_number, envir = envProj)
                                    assign("locus_set", locus_set, envir = envProj)
                                    assign("locus_set_name", locus_set_name, envir = envProj)
                                    assign("LinkMode", tclvalue(LinkModeVar), envir = envProj)

                                    tkdestroy(tf)
                                  })



  label_Linkage.locus1 <- tklabel(frame_Linkage.label,text="Locus 1",width="12")
  label_Linkage.locus2 <- tklabel(frame_Linkage.label,text="Locus 2",width="12")
  label_Linkage.RR <- tklabel(frame_Linkage.label,text="Recombination rate",width="17")


  tkpack(outframe_Linkage.setting, anchor="w", padx=40, pady="20")
  tkpack(frame_Linkage.radiobutton, anchor="w", padx=20, pady="15")

  if (tclvalue(LinkModeVar) == "1"){
    tkpack(frame_Linkage.add, anchor="w", padx=40, pady=5)
    tkpack(frame_Linkage.label, anchor="w", padx=40, pady=5)
    tkpack(frame_Linkage.select, padx=40, anchor="w")
    if (linkage_number > 0){
      for (i in 1:linkage_number){
        tkpack.forget(frame_Linkage.select.input[[i]])
      }
    }
    linkage_number.default <<- length(RR.default)
    default.sort <- sapply(default_locus_set[[locus_set_name]][[1]], function(x){return(which(x == locus_set))})
    for (i in 1:linkage_number.default){
      tkpack(frame_Linkage.select.input[[i]], anchor="w", pady="5")
      tkconfigure(combobox_Linkage.locus1[[i]], state="disable")
      tkconfigure(combobox_Linkage.locus2[[i]], state="disable")
      tkconfigure(entry_Linkage.RR[[i]], state="disable",cursor="arrow")
      tkconfigure(button_Linkage.vanish[[i]], state="disable",cursor="arrow")
      tclvalue(locus1.var[[i]]) <- var_locus[default.sort[linkage.loci.default[i,1]]]
      tclvalue(locus2.var[[i]]) <- var_locus[default.sort[linkage.loci.default[i,2]]]
      tclvalue(RR.var[[i]]) <- RR.default[i]
    }
    tkconfigure(button_Linkage.add, state="disable",cursor="arrow")

  }else if(tclvalue(LinkModeVar) == "2"){
    tkpack(frame_Linkage.add, anchor="w", padx=40, pady=5)
    tkpack(frame_Linkage.label, anchor="w", padx=40, pady=5)
    tkpack(frame_Linkage.select, padx=40, anchor="w")
    for (i in 1:pair.max){

      tkconfigure(combobox_Linkage.locus1[[i]], state="normal")
      tkconfigure(combobox_Linkage.locus2[[i]], state="normal")
      tkconfigure(entry_Linkage.RR[[i]], state="normal",cursor="arrow")
      tkconfigure(button_Linkage.vanish[[i]], state="normal",cursor="arrow")
    }
    tkpack(frame_Linkage.select.input[[1]], anchor="w", pady="5")
    tkconfigure(button_Linkage.add, state="normal",cursor="hand2")
  }else{
    tkpack.forget(frame_Linkage.add)
    tkpack.forget(frame_Linkage.label)
    tkpack.forget(frame_Linkage.select)
  }
  tkpack(frame_Linkage.save, fill="both", expand="1", anchor="w", padx=10, pady=10)

  tkpack(label_Linkage.locus1, side="left", padx=5)
  tkpack(label_Linkage.locus2, side="left", padx=5)
  tkpack(label_Linkage.RR, anchor="w", side="left", padx=5)

  for (i in 1:pair.max){
    tkpack(combobox_Linkage.locus1[[i]], padx=5, side="left")
    tkpack(combobox_Linkage.locus2[[i]], padx=5, side="left")
    tkpack(entry_Linkage.RR[[i]], padx=5, side="left")
    tkpack(button_Linkage.vanish[[i]], padx=5, side="left")
  }

  if (locus_set_name != "-"){
    tkpack(radiobutton_off, anchor="w", padx=15)
  }
  tkpack(radiobutton_on, anchor="w", padx=15)
  tkpack(radiobutton_none, anchor="w", padx=15)

  tkpack(button_Linkage.add, padx=10, side="left")
  tkpack(frame_Linkage.save)
  tkpack(button_Linkage.save, anchor="e", padx=15)
}

