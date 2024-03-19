makeHypo <- function(envProj, envGUI, modeStateVar){
  tabHypo <- get("tabHypo", pos = envGUI)
  tabResults <- get("tabResults", pos = envGUI)
  finHypo <- get("finHypo", pos = envProj)
  frame_Hypo <- tkframe(tabHypo)
  frame_Hypo.apply <- tkframe(tabHypo)

  if(finHypo == FALSE){
    labelframe_Hypo <- list()
    frame_Hypo.button <- list()
    frame_Hypo.label <- list()
    for (j in 1:2){
      labelframe_Hypo[[j]] <- tk2labelframe(tabHypo,text=paste("Hypothesis ", j, " (H", j, ")", collapse="", sep=""),labelanchor="nw",relief="groove",borderwidth="2")
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

    checkFile(envProj, modeStateVar)
    linkState <- get("linkState", pos = envProj)
    dropState <- get("dropState", pos = envProj)
    freqState <- get("freqState", pos = envProj)
    setting_menu <- get("setting_menu", pos = envGUI)

    if(linkState=="normal"){

      linkStateVar <- tclVar(linkState)
      assign("linkState", tclVar(linkStateVar), envir = envProj)
      dropStateVar <- tclVar(dropState)
      assign("dropState", tclVar(dropStateVar), envir = envProj)
      freqStateVar <- tclVar(freqState)
      assign("freqState", tclVar(freqStateVar), envir = envProj)
      tkdelete(setting_menu, 0, 2)

      tkadd(setting_menu, "command", state = tclvalue(linkStateVar), label = "Linkage", command = function() windowLinkSet(envProj))
      tkadd(setting_menu, "command", state = tclvalue(dropStateVar), label = "Drop-out", command = function() windowDropSet(envProj))
      tkadd(setting_menu, "command", state = tclvalue(freqStateVar), label = "Allele frequency", command = function() windowFreqSet(envProj))
    }

    ped.max <- get("ped.max", pos = envProj)
    DNAtype_matrix <- get("DNAtype_matrix", pos = envProj)
    labelframe_Hypo <- get("labelframe_Hypo", pos = envGUI)
    frame_Hypo.button <- get("frame_Hypo.button", pos = envGUI)
    frame_Hypo.label <- get("frame_Hypo.label", pos = envGUI)
    frame_Hypo.next <- get("frame_Hypo.next", pos = envGUI)
    frame_Hypo.select.out <- get("frame_Hypo.select.out", pos = envGUI)
    canvas_Hypo <- get("canvas_Hypo", pos = envGUI)
    scr_Hypo <- get("scr_Hypo", pos = envGUI)
    frame_Hypo.select <- get("frame_Hypo.select", pos = envGUI)

    if (modeStateVar == "1"){
      var_ref <-  names(DNAtype_matrix)[1:((length(DNAtype_matrix) - 1) / 2) * 2]
      var_Hypo <- c(var_ref, paste("UK",1:(ped.max - length(var_ref)),sep=""))
    }else{
      var_ref <- ""
      var_Hypo <- rep("", ped.max)
    }

    var_sex <- c("Male", "Female") #sex
    var_parent <- list(var_ref, var_ref)

    sex.var <- list(lapply(1:length(var_Hypo),function(i){tclVar("")}), lapply(1:length(var_Hypo),function(i){tclVar("")}))
    pat.var <- list(lapply(1:length(var_Hypo),function(i){tclVar("")}), lapply(1:length(var_Hypo),function(i){tclVar("")}))
    mat.var <- list(lapply(1:length(var_Hypo),function(i){tclVar("")}), lapply(1:length(var_Hypo),function(i){tclVar("")}))
    founder.var <- list(lapply(1:length(var_Hypo),function(i){tclVar("0")}), lapply(1:length(var_Hypo),function(i){tclVar("0")}))
    ped.number <- rep(length(var_ref), 2)
    N <- tclVar("")
    ref.num <- tclVar("")

    if (modeStateVar == "2"){

      label_ref.num <- tklabel(frame_Hypo.apply, text="Number of known profiles", anchor="w")
      text_ref.num <- tkentry(frame_Hypo.apply, textvariable=ref.num, width=8, background="#ffffff")
      button_ref.num <- tkbutton(frame_Hypo.apply, text="Apply", cursor="hand2", width=8,
                                 command = function(){
                                   if (is.na(as.integer(tclvalue(ref.num)))){
                                     tkmessageBox(message = paste("Please enter an integer of 2-", ped.max, " to Number of known profiles.", sep="", collapse=""), icon="error")
                                   }else if ((as.integer((tclvalue(ref.num))) < 2 || as.integer(tclvalue(ref.num)) > ped.max) ){
                                     tkmessageBox(message = paste("Please enter an integer of 2-", ped.max, " to Number of known profiles.", sep="", collapse=""), icon="error")
                                   }else{
                                     for (j in 1:2){
                                       for (i in 1:ped.number[j]){
                                         tkpack.forget(frame_Hypo.select.input[[j]][[i]])
                                         tclvalue(sex.var[[j]][[i]]) <<- ""
                                         tclvalue(pat.var[[j]][[i]]) <<- ""
                                         tclvalue(mat.var[[j]][[i]]) <<- ""
                                         tclvalue(founder.var[[j]][[i]]) <<- 0
                                         tkconfigure(combobox_Hypo.pat[[j]][[i]], state="readonly")
                                         tkconfigure(combobox_Hypo.mat[[j]][[i]], state="readonly")
                                       }
                                     }

                                     var_ref <<- paste("ID", 1:as.integer((tclvalue(ref.num))), sep="")
                                     var_Hypo <<- c(paste("ID", 1:as.integer((tclvalue(ref.num))), sep=""), paste("UK",1:(ped.max - length(var_ref)),sep=""))
                                     var_parent <<- list(var_ref, var_ref)
                                     tkpack(labelframe_Hypo[[1]], anchor="n", padx="20", pady="10", side="left")
                                     tkpack(labelframe_Hypo[[2]], anchor="n", padx="20", pady="10")
                                     tkpack(frame_Hypo.next, anchor="e", padx="10", pady="10")
                                     for (j in 1:2){
                                       for (i in 1:ped.max){
                                         tkconfigure(combobox_Hypo.sex[[j]][[i]], textvariable=sex.var[[j]][[i]])
                                         tkconfigure(label_Hypo.listname[[j]][[i]], text=as.list(var_Hypo)[[i]])
                                         tkconfigure(combobox_Hypo.pat[[j]][[i]], values=var_parent[[j]])
                                         tkconfigure(combobox_Hypo.mat[[j]][[i]], values=var_parent[[j]])
                                       }
                                       for (i in 1:length(var_ref)){
                                         tkconfigure(combobox_Hypo.sex[[j]][[i]], textvariable=sex.var[[1]][[i]])
                                       }
                                       for (i in 1:as.integer((tclvalue(ref.num)))){
                                         tkpack(frame_Hypo.select.input[[j]][[i]],anchor="w",padx="10",pady="3")
                                       }
                                       ped.number[j] <<- as.integer((tclvalue(ref.num)))
                                     }
                                     tkconfigure(canvas_Hypo[[1]], scrollregion=c(0,0,1000,32 * length(var_ref)))
                                     tkconfigure(canvas_Hypo[[2]], scrollregion=c(0,0,1000,32 * length(var_ref)))
                                   }
                                 })
    }
    checkbutton_Hypo.change <- vector("list", length = 2)
    checkbutton_Hypo.change[[1]] <- function(){
      for (i in 1:ped.max){
        if (tclvalue(founder.var[[1]][[i]]) == 1){
          tkconfigure(combobox_Hypo.pat[[1]][[i]], state="disable",cursor="arrow")
          tkconfigure(combobox_Hypo.mat[[1]][[i]], state="disable",cursor="arrow")
          tclvalue(pat.var[[1]][[i]]) <- ""
          tclvalue(mat.var[[1]][[i]]) <- ""
        }else{
          tkconfigure(combobox_Hypo.pat[[1]][[i]], state="readonly")
          tkconfigure(combobox_Hypo.mat[[1]][[i]], state="readonly")
        }
      }
    }
    checkbutton_Hypo.change[[2]] <- function(){
      for (i in 1:ped.max){
        if (tclvalue(founder.var[[2]][[i]]) == 1){
          tkconfigure(combobox_Hypo.pat[[2]][[i]], state="disable",cursor="arrow")
          tkconfigure(combobox_Hypo.mat[[2]][[i]], state="disable",cursor="arrow")
          tclvalue(pat.var[[2]][[i]]) <- ""
          tclvalue(mat.var[[2]][[i]]) <- ""
        }else{
          tkconfigure(combobox_Hypo.pat[[2]][[i]], state="readonly")
          tkconfigure(combobox_Hypo.mat[[2]][[i]], state="readonly")
        }
      }
    }

    frame_Hypo.select.input <- label_Hypo.listname <- combobox_Hypo.sex <- combobox_Hypo.pat <- combobox_Hypo.mat <- checkbutton_Hypo.founder <- vector("list", length = 2)

    for (j in 1:2){
      for (i in 1:ped.max){
        frame_Hypo.select.input[[j]][[i]] <- tkframe(frame_Hypo.select[[j]])
        label_Hypo.listname[[j]][[i]] <- tklabel(frame_Hypo.select.input[[j]][[i]], text=as.list(var_Hypo)[[i]], width=12)
        combobox_Hypo.sex[[j]][[i]] <- ttkcombobox(frame_Hypo.select.input[[j]][[i]], values=var_sex, textvariable=sex.var[[j]][[i]], cursor="hand2", width=7, state="readonly")
        combobox_Hypo.pat[[j]][[i]] <- ttkcombobox(frame_Hypo.select.input[[j]][[i]], values=var_parent[[j]], textvariable=pat.var[[j]][[i]], cursor="hand2", width=10, state="readonly")
        combobox_Hypo.mat[[j]][[i]] <- ttkcombobox(frame_Hypo.select.input[[j]][[i]], values=var_parent[[j]], textvariable=mat.var[[j]][[i]], cursor="hand2", width=10, state="readonly")
        checkbutton_Hypo.founder[[j]][[i]] <- tkcheckbutton(frame_Hypo.select.input[[j]][[i]], variable=founder.var[[j]][[i]], width=5, command=checkbutton_Hypo.change[[j]],cursor="hand2")
      }
    }
    if (modeStateVar == "1"){
      for (i in 1:length(var_ref)){
        tkconfigure(combobox_Hypo.sex[[j]][[i]], textvariable=sex.var[[1]][[i]])
      }
    }

    Hypo.add <- function(j){
      if (length(var_ref) <= ped.number[j] && ped.number[j] < ped.max){
        ped.number[j] <<- ped.number[j] + 1
        var_parent[[j]] <<- c(var_parent[[j]], var_Hypo[ped.number[j]])
        for (i in 1:ped.max){
          tkconfigure(combobox_Hypo.pat[[j]][[i]], values=var_parent[[j]])
          tkconfigure(combobox_Hypo.mat[[j]][[i]], values=var_parent[[j]])
        }
        tkpack(frame_Hypo.select.input[[j]][[ped.number[j]]], anchor="w", padx="10", pady="3")
        tkconfigure(canvas_Hypo[[j]], scrollregion=c(0,0,0,32 * ped.number[j]))
      }
    }

    Hypo.vanish <- function(j){
      if (length(var_ref) < ped.number[j] && ped.number[j] <= ped.max){
        tkpack.forget(frame_Hypo.select.input[[j]][[ped.number[j]]])
        tclvalue(sex.var[[j]][[ped.number[j]]]) <<- ""
        tclvalue(pat.var[[j]][[ped.number[j]]]) <<- ""
        tclvalue(mat.var[[j]][[ped.number[j]]]) <<- ""
        tclvalue(founder.var[[j]][[ped.number[j]]]) <<- 0
        tkconfigure(combobox_Hypo.pat[[j]][[ped.number[j]]], state="readonly")
        tkconfigure(combobox_Hypo.mat[[j]][[ped.number[j]]], state="readonly")
        var_parent[[j]] <<- var_parent[[j]][-ped.number[j]]
        for (i in 1:ped.max){
          tkconfigure(combobox_Hypo.pat[[j]][[i]], values=var_parent[[j]])
          tkconfigure(combobox_Hypo.mat[[j]][[i]], values=var_parent[[j]])
        }
        ped.number[j] <<- ped.number[j] - 1
        tkconfigure(canvas_Hypo[[j]], scrollregion=c(0,0,0,32 * ped.number[j]))
      }
    }

    sex.var_to_val <- function(x){
      if (tclvalue(x) == "Male"){
        return(1)
      }else if (tclvalue(x) == "Female"){
        return(2)
      }else{
        return(-1)
      }
    }

    parent.var_to_val <- function(x, var, founder, ped.number){
      ID <- which(c(var_Hypo, "") == tclvalue(var[[x]])) * (1 - as.integer(tclvalue(founder[[x]])))
      if (ID > ped.number){
        return(-1)
      }else{
        return(ID)
      }
    }

    Hypo.pedigree.check <- function(j){
      if (j == 1){
        sex <- sapply(sex.var[[1]][1:ped.number[1]],sex.var_to_val)
      }else{
        sex <- sapply(sex.var[[1]][1:length(var_ref)],sex.var_to_val)
        if (ped.number[2] != length(var_ref)){
          sex <- c(sex, sapply(sex.var[[2]][(length(var_ref) + 1) : ped.number[2]],sex.var_to_val))
        }
      }
      pat <- sapply(1:ped.number[j], parent.var_to_val, pat.var[[j]], founder.var[[j]], ped.number[j])
      mat <- sapply(1:ped.number[j], parent.var_to_val, mat.var[[j]], founder.var[[j]], ped.number[j])
      affected <- c(rep(1,length(var_ref)),rep(0, (ped.number[j] - length(var_ref))))
      DNAtype <- c(1:length(var_ref),rep(0, (ped.number[j] - length(var_ref))))
      df <- data.frame(id = 1:length(var_parent[[j]]), "names" = var_parent[[j]], "sex" = sex, "pat" = pat, "mat" = mat, "affected" = affected, "DNAtype" = DNAtype)
      if (is.element(-1, df[["sex"]])){
        tkmessageBox(message = paste("Hypothesis", j, ": Error: incorrect setting of Sex", collapse="", sep=""), icon="error")
        stop()
      }else if (is.element(-1, df$pat) || is.element(-1, df$mat)){
        tkmessageBox(message = paste("Hypothesis", j, ": Error: incorrect setting of Father or Mother", collapse="", sep=""), icon="error")
        stop()
      }else{
        err <- try(pedigree(id = df$id, sex= df$sex, dadid = df$pat, momid = df$mat), silent = T)
        err.fc <- try(familycheck(famid=rep(1,nrow(df)), id = df$id, father.id = df$pat, mother.id = df$mat), silent = T)
        if ( !is.na(match(class(err), "try-error")) || !is.na(match(class(err.fc), "try-error"))) {
          tkmessageBox(message = paste("Hypothesis", j, ": Error: incorrect setting of pedigree ", collapse="", sep=""), icon="error")
          stop()
        }else{
          return(df)
        }
      }
    }

    Hypo.pedigree <- function(df, x){
      if (is.data.frame(df)){
        unrelated <- which(df$pat == 0)[is.element(which(df$pat == 0), c(df$pat, df$mat)) == F]
        if (length(unrelated) != nrow(df)){
          if (length(unrelated) != 0){
            pedobj <- pedigree(id = df[-unrelated,]$id, sex= df[-unrelated,]$sex, dadid = df[-unrelated,]$pat, momid = df[-unrelated,]$mat, affected = df[-unrelated,]$affected)
            plot(pedobj, id = df[-unrelated,]$names, main = paste("Hypothesis", x, sep = " "))
          }else{
            pedobj <- pedigree(id = df$id, sex= df$sex, dadid = df$pat, momid = df$mat, affected = df$affected)
            plot(pedobj, id = df$names, main = paste("Hypothesis", x, sep = " "))
          }
        }
        if (length(unrelated) != 0){
          tkmessageBox(title = paste("Hypothesis", x, sep = " "),
                       message = paste("Unrelated", paste(df[unrelated,]$names, collapse=", "), sep = " : "),
                       icon = "info")
        }
      }
    }

    button_Hypo.add <- button_Hypo.vanish <- button_Hypo.pedigree <- vector("list", length = 2)
    button_Hypo.add[[1]] <- tkbutton(frame_Hypo.button[[1]], text="Add", cursor="hand2", width=9, command = function(){Hypo.add(1)} )
    button_Hypo.add[[2]] <- tkbutton(frame_Hypo.button[[2]], text="Add", cursor="hand2", width=9, command = function(){Hypo.add(2)} )
    button_Hypo.vanish[[1]] <- tkbutton(frame_Hypo.button[[1]], text="Delete", cursor="hand2", width=9, command = function(){Hypo.vanish(1)} )
    button_Hypo.vanish[[2]] <- tkbutton(frame_Hypo.button[[2]], text="Delete", cursor="hand2", width=9, command = function(){Hypo.vanish(2)} )
    button_Hypo.pedigree[[1]] <- tkbutton(frame_Hypo.button[[1]], text="View pedigree tree",cursor="hand2", command = function(){Hypo.pedigree(Hypo.pedigree.check(1), 1)} )
    button_Hypo.pedigree[[2]] <- tkbutton(frame_Hypo.button[[2]], text="View pedigree tree",cursor="hand2", command = function(){Hypo.pedigree(Hypo.pedigree.check(2), 2)} )

    label_Hypo.name <- label_Hypo.sex <- label_Hypo.pat <- label_Hypo.mat <- label_Hypo.founder <- vector("list", length = 2)

    for (j in 1:2){
      label_Hypo.name[[j]] <- tklabel(frame_Hypo.label[[j]],text="Name",  width=12)
      label_Hypo.sex[[j]] <- tklabel(frame_Hypo.label[[j]],text="Sex",  width=8)
      label_Hypo.pat[[j]] <- tklabel(frame_Hypo.label[[j]],text="Father",  width=12)
      label_Hypo.mat[[j]] <- tklabel(frame_Hypo.label[[j]],text="Mother",  width=12)
      label_Hypo.founder[[j]] <- tklabel(frame_Hypo.label[[j]],text="founder", width=7)
    }
    if (modeStateVar == "2"){
      N <- tclVar("")
      label_N <- tklabel(frame_Hypo.next, text="Number of simulations", width=18)
      text_N <- tkentry(frame_Hypo.next, textvariable=N, width=8,  background="#ffffff")
    }
    button_Hypo.next <- tkbutton(frame_Hypo.next, text="Calculate", cursor="hand2", width=10,
                                 command = function() {
                                   assign("df_Hypo1", Hypo.pedigree.check(1), envir = envProj)
                                   assign("df_Hypo2", Hypo.pedigree.check(2), envir = envProj)
                                   calcLR(envProj, envGUI)
                                 })
    button_Hypo.next.simu <- tkbutton(frame_Hypo.next, text="Calculate", cursor="hand2", width=10,
                                      command = function() {
                                        assign("df_Hypo1", Hypo.pedigree.check(1), envir = envProj)
                                        assign("df_Hypo2", Hypo.pedigree.check(2), envir = envProj)
                                        calcSimu(envProj, envGUI, tclvalue(N), modeStateVar)
                                      })
    if (modeStateVar == "1"){
      tkpack(labelframe_Hypo[[1]], anchor="n", padx=20, pady=10, side="left")
      tkpack(labelframe_Hypo[[2]], anchor="n", padx=20, pady=10)
      tkpack(frame_Hypo.next, anchor="e", padx=10, pady=10)
    }else{
      tkpack(frame_Hypo.apply, anchor="w", padx=20, pady=20)
      tkpack(label_ref.num, anchor="w", padx=5, side="left")
      tkpack(text_ref.num, padx=5, side="left")
      tkpack(button_ref.num, padx=10)
    }

    for (j in 1:2){
      tkpack(frame_Hypo.button[[j]],anchor="w",padx="10",pady="10")
      tkpack(frame_Hypo.label[[j]],anchor="w",padx="10",pady="10")
      tkpack(frame_Hypo.select.out[[j]], anchor="w")
      if (modeStateVar == "1"){
        for (i in 1:length(var_ref)){
          tkpack(frame_Hypo.select.input[[j]][[i]],anchor="w",padx="10",pady="3")
        }
      }

      tkpack(label_Hypo.name[[j]], side="left", padx=5)
      tkpack(label_Hypo.sex[[j]], side="left", padx=5)
      tkpack(label_Hypo.pat[[j]], side="left", padx=5)
      tkpack(label_Hypo.mat[[j]], side="left", padx=5)
      tkpack(label_Hypo.founder[[j]], padx=5)

      tkpack(button_Hypo.add[[j]], padx=10, side="left")
      tkpack(button_Hypo.vanish[[j]], padx=10, side="left")
      tkpack(button_Hypo.pedigree[[j]], padx=10)

      for (i in 1:ped.max){
        tkpack(label_Hypo.listname[[j]][[i]], padx=5, side="left")
        tkpack(combobox_Hypo.sex[[j]][[i]], padx=5, side="left")
        tkpack(combobox_Hypo.pat[[j]][[i]], padx=5, side="left")
        tkpack(combobox_Hypo.mat[[j]][[i]], padx=5, side="left")
        tkpack(checkbutton_Hypo.founder[[j]][[i]], anchor="w", padx=5)
      }
    }
    tkpack(scr_Hypo[[1]])
    tkconfigure(scr_Hypo[[1]], command=function(...)tkyview(canvas_Hypo[[1]],...))
    tkpack.configure(scr_Hypo[[1]], side="right", expand="1", fill="y")
    tkpack(canvas_Hypo[[1]], side="left")
    tkconfigure(canvas_Hypo[[1]], width=505, height=292, scrollregion=c(0,0,1000,32 * length(var_ref)), yscrollcommand=function(...)tkset(scr_Hypo[[1]],...))
    tkcreate(canvas_Hypo[[1]], "window", 0, 0, anchor="nw", window=frame_Hypo.select[[1]])
    tkpack(scr_Hypo[[2]])
    tkconfigure(scr_Hypo[[2]], command=function(...)tkyview(canvas_Hypo[[2]],...))
    tkpack.configure(scr_Hypo[[2]], side="right", expand="1", fill="y")
    tkpack(canvas_Hypo[[2]], side="left")
    tkconfigure(canvas_Hypo[[2]], width=505, height=292, scrollregion=c(0,0,1000,32 * length(var_ref)), yscrollcommand=function(...)tkset(scr_Hypo[[2]],...))
    tkcreate(canvas_Hypo[[2]], "window", 0, 0, anchor="nw", window=frame_Hypo.select[[2]])

    if (modeStateVar == "2"){
      tkpack(label_N, padx="5", pady="5", anchor="e", side="left")
      tkpack(text_N, padx="5", pady="5", anchor="e", side="left")
      tkpack(button_Hypo.next.simu, anchor="e", padx="15")
    }else{
      tkpack(button_Hypo.next, anchor="e", padx="15")
    }
  }
  assign("frame_Hypo.apply", frame_Hypo.apply, envir = envGUI)

  ref.numVar <- tclVar(ref.num)
  assign("ref.num", tclvalue(ref.numVar), envir = envProj)
  tabs <- get("tabs", pos = envGUI)
  tk2notetab.select(tabs, "Hypothesis")
  assign("finHypo", TRUE, envir = envProj)

}
