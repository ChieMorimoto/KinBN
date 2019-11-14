KinBN <- function(ped.max = 25, range = 0){
  version <- "ver 1.1.0"
  
  if(require("bnlearn")){
    print("bnlearn is loaded correctly")
  }
  else{
    print("trying to install bnlearn...")
    install.packages("bnlearn")
    if(require("bnlearn")){
      print("bnlearn installed and loaded")
    }
    else{
      stop("could not install bnlearn")
    }
  }

  if(require("BiocManager")){
    print("BiocManager is loaded correctly")
  }
  else{
    print("trying to install BiocManager...")
    install.packages("BiocManager")
    if(require("BiocManager")){
      print("BiocManager installed and loaded")
    }
    else{
      stop("could not install BiocManager")
    }
  }

  if(requireNamespace("BiocGenerics")){
    print("BiocGenerics is loaded correctly")
  }
  else{
    print("trying to install BiocGenerics...")
    BiocManager::install(c("BiocGenerics"),update=F)
    if(requireNamespace("BiocGenerics")){
      print("BiocGenerics installed and loaded")
    }
    else{
      stop("could not install BiocGenerics")
    }
  }

  if(requireNamespace("RBGL")){
    print("RBGL is loaded correctly")
  }
  else{
    print("trying to install RBGL...")
    BiocManager::install(c("RBGL"),update=F)
    if(requireNamespace("RBGL")){
      print("RBGL installed and loaded")
    }
    else{
      stop("could not install RBGL")
    }
  }

  if(requireNamespace("graph")){
    print("graph is loaded correctly")
  }
  else{
    print("trying to install graph...")
    source("https://bioconductor.org/biocLite.R")
    biocLite(c("graph"),suppressUpdates=TRUE)
    if(requireNamespace("graph")){
      print("graph installed and loaded")
    }
    else{
      stop("could not install graph")
    }
  }

  if(require("gRbase")){
    print("gRbase is loaded correctly")
  }
  else{
    print("trying to install gRbase...")
   install.packages("gRbase")
    if(require("gRbase")){
      print("gRbase installed and loaded")
    }
    else{
      stop("could not install gRbase")
    }
  }

  if(require("gRain")){
    print("gRain is loaded correctly")
  }
  else{
    print("trying to install gRain...")
    install.packages("gRain")
    if(require("gRain")){
      print("gRain installed and loaded")
    }
    else{
      stop("could not install gRain")
    }
  }

  if(require("tcltk")){
    print("tcltk is loaded correctly")
  }
  else{
    print("trying to install tcltk...")
    install.packages("tcltk")
    if(require("tcltk")){
      print("tcltk installed and loaded")
    }
    else{
      stop("could not install tcltk")
    }
  }

  if(require("tcltk2")){
    print("tcltk2 is loaded correctly")
  }
  else{
    print("trying to install tcltk2...")
    install.packages("tcltk2")
    if(require("tcltk2")){
      print("tcltk2 installed and loaded")
    }
    else{
      stop("could not install tcltk2")
    }
  }

  if(require("tkRplotR")){
    print("tkRplotR is loaded correctly")
  }
  else{
    print("trying to install tkRplotR...")
    install.packages("tkRplotR")
    if(require("tkRplotR")){
      print("tkRplotR installed and loaded")
    }
    else{
      stop("could not install tkRplotR")
    }
  }

  if(require("kinship2")){
    print("kinship2 is loaded correctly")
  }
  else{
    print("trying to install kinship2...")
    install.packages("kinship2")
    if(require("kinship2")){
      print("kinship2 installed and loaded")
    }
    else{
      stop("could not install kinship2")
    }
  }

  OpenFile <- function(fp,check,top){
    if (!already_pack.Linkage_Hypo || tclvalue(tkmessageBox(message="Setting of Linkage and Hypothesis will be deleted. Do you want to continue?", type="okcancel", parent = window, icon="warning")) == "ok"){
      if (already_pack.Linkage_Hypo){
        tkdestroy(frame_Linkage)
        tkdestroy(frame_Hypo)
        tkdestroy(frame_LR_result)
        tkdestroy(frame_LR_simu)
        tkdestroy(frame_Simu)
        already_pack.Linkage_Hypo <<- already_pack.apply <<- already_pack.LR <<- already_pack.simu <<- F
      }
      FileName <- tclvalue(tkgetOpenFile(parent = top, multiple = "false", filetypes = "{{CSV Files} {.csv}}"))
      if (nchar(FileName)){
        tclvalue(fp) <- FileName
        tclvalue(check) <- strsplit(FileName,"/")[[1]][length(strsplit(FileName,"/")[[1]])]
      }
    }
  }

  default_locus_set <- list() 

  ##GlobalFiler
  GF.set <- c("D3S1358","vWA","D16S539","CSF1PO","TPOX","D8S1179","D21S11","D18S51","D2S441","D19S433","TH01","FGA","D22S1045","D5S818","D13S317","D7S820","SE33","D10S1248","D1S1656","D12S391","D2S1338")
  GF.linkage.loci <- rbind(c(9,5),c(4,14),c(20,2))
  GF.RR <- c(0.4721,0.2522,0.1172)
  default_locus_set[[1]] <- list(GF.set,GF.linkage.loci,GF.RR)
  names(default_locus_set)[1] <- "GlobalFiler"
  ##Identifiler
  ID.set <- c("D8S1179", "D21S11", "D7S820", "CSF1PO", "D3S1358", "TH01", "D13S317", "D16S539", "D2S1338", "D19S433", "vWA", "TPOX", "D18S51", "D5S818", "FGA")
  ID.linkage.loci <- rbind(c(4,14))
  ID.RR <- c(0.2522)
  default_locus_set[[2]] <- list(ID.set,ID.linkage.loci,ID.RR)
  names(default_locus_set)[2] <- "Identifiler"
  ##PowerPlex Fusion
  PPF.set <- c("D3S1358","D1S1656","D2S441","D10S1248","D13S317","Penta E","D16S539","D18S51","D2S1338","CSF1PO","Penta D","TH01","vWA","D21S11","D7S820","D5S818","TPOX","D8S1179","D12S391","D19S433","FGA","D22S1045")
  PPF.linkage.loci <- rbind(c(20,13),c(16,10),c(14,11),c(3,17))
  PPF.RR <- c(0.1172,0.2522,0.3568,0.4721)
  default_locus_set[[3]] <- list(PPF.set,PPF.linkage.loci,PPF.RR)
  names(default_locus_set)[3] <- "PowerPlex Fusion"
  ##PowerPlex 16
  PP16.set <- c("D3S1358","TH01","D21S11","D18S51","Penta E","D5S818","D13S317","D7S820","D16S539","CSF1PO","Penta D","vWA","D8S1179","TPOX","FGA")
  PP16.linkage.loci <- rbind(c(6,10),c(3,11))
  PP16.RR <- c(0.2522,0.3568)
  default_locus_set[[4]] <- list(PP16.set,PP16.linkage.loci,PP16.RR)
  names(default_locus_set)[4] <- "PowerPlex 16"
  ##NGM SElect
  NGMSE.set <- c("D10S1248","vWA","D16S539","D2S1338","D8S1179","D21S11","D18S51","D22S1045","D19S433","TH01","FGA","D2S441","D3S1358","D1S1656","D12S391","SE33")
  NGMSE.linkage.loci <- rbind(c(15,2))
  NGMSE.RR <- c(0.1172)
  default_locus_set[[5]] <- list(NGMSE.set,NGMSE.linkage.loci,NGMSE.RR)
  names(default_locus_set)[5] <- "NGM SElect"

  checkFile <- function(cf){
    if (tclvalue(cf$filepath_frequency) == ""){
      tkmessageBox(message = "Please load allele frequency file.", parent = window, icon="error")
    }
    else if(tclvalue(cf$filepath_mutation) == ""){
      tkmessageBox(message = "Please load mutation rates file.", parent = window, icon="error")
    }
    else if(cf$mode == "LR" && tclvalue(cf$filepath_DNAtype) == ""){
      tkmessageBox(message = "Please load profiles file.", parent = window, icon="error")
    }
    else{
      frequency_matrix <- read.csv(tclvalue(cf$filepath_frequency))
      locus_set <- names(frequency_matrix)[-1]
      mutation_matrix <- read.csv(tclvalue(cf$filepath_mutation), stringsAsFactors=F )
      locus_set.mutation <- mutation_matrix[,1]
      if (cf$mode == "LR"){
        DNAtype_matrix <- read.csv(tclvalue(cf$filepath_DNAtype), stringsAsFactors=F )
        locus_set.DNAtype <- DNAtype_matrix[,1]
      }
      else{
        DNAtype_matrix <- ""
      }

      if (ncol(mutation_matrix) != 11){
        tkmessageBox(message = "Error: incorrect file format of Mutation rate", parent = window, icon="error")
      }
      else if (!identical(locus_set, locus_set.mutation)){
        tkmessageBox(message = "Error: incorrect locus set", parent = window, icon=="error")
      }
      else if (cf$mode == "LR" && prod(sapply(locus_set, function(x){return(sum(locus_set.DNAtype == x) == 1)})) == 0 ){
        tkmessageBox(message = "Error: incorrect locus set", parent = window, icon="error")
      }
      else{
        is.default <- sapply(default_locus_set, function(x){if (setequal(x[[1]], locus_set) && length(x[[1]]) == length(locus_set) ) return(1) else return(0)})
        if (length(which(is.default == 1)) > 0){
          locus_set_name <- names(default_locus_set)[which(is.default == 1)]
        }
        else{
          locus_set_name <- "-"
        }

        if (cf$mode == "LR"){
          DNAtype.sort <- sapply(locus_set, function(x){return(which(x == locus_set.DNAtype))})
          DNAtype_matrix <- DNAtype_matrix[DNAtype.sort, ]
        }

        if (cf$mode == "simu" || length(locus_set) == length(locus_set.DNAtype) || tclvalue(tkmessageBox(message="Loci without allele frequencies and mutation rates will be excluded from calculation.\nDo you want to continue?", type="okcancel", parent = window, icon="warning")) == "ok"){
          temp <- list(frequency_matrix, mutation_matrix, DNAtype_matrix, locus_set_name, locus_set)
          names(temp) <- c("frequency_matrix", "mutation_matrix", "DNAtype_matrix", "locus_set_name", "locus_set")
          cf.result <- c(cf, temp)
          return(cf.result)
        }
      }
    }
  }

  makeFiles <- function(mode){
    filepath_frequency <- tclVar("")
    filepath_DNAtype <- tclVar("")
    filepath_mutation <- tclVar("")
    OK_frequency <- tclVar("")
    OK_DNAtype <- tclVar("")
    OK_mutation <- tclVar("")
    ref.num <- tclVar("")

    outframe.file <- tkframe(frame_Files)
    frame_file.frequency <- tkframe(outframe.file)
    frame_file.mutation <- tkframe(outframe.file)
    if (mode == "LR") frame_file.DNAtype <- tkframe(outframe.file)
    frame_Files.next <- tkframe(frame_Files)
  
    label_frequency <- tklabel(frame_file.frequency, text="Allele frequencies", font="our_font", width=15, anchor="w")
    check_frequency <- tklabel(frame_file.frequency, textvariable=OK_frequency, font="our_font", width=30, relief="groove", borderwidth="3", anchor="w")
    button_frequency <- tkbutton(frame_file.frequency, text="Browse", font="our_font", width=10, command=function() OpenFile(filepath_frequency, OK_frequency, window) )
    label_mutation <- tklabel(frame_file.mutation, text="Mutation rates", font="our_font", width=15, anchor="w")
    check_mutation <- tklabel(frame_file.mutation, textvariable=OK_mutation, font="our_font", width=30, relief="groove", borderwidth="3", anchor="w")
    button_mutation <- tkbutton(frame_file.mutation, text="Browse", font="our_font", width=10, command=function() OpenFile(filepath_mutation, OK_mutation, window) )
    if (mode == "LR"){
      label_DNAtype <- tklabel(frame_file.DNAtype, text="Profiles", font="our_font",width=15,anchor="w")
      check_DNAtype <- tklabel(frame_file.DNAtype, textvariable=OK_DNAtype, font="our_font", width=30, relief="groove", borderwidth="3", anchor="w")
      button_DNAtype <- tkbutton(frame_file.DNAtype, text="Browse", font="our_font", width=10, command=function() OpenFile(filepath_DNAtype, OK_DNAtype, window) )
    }
    button_Files.next <- tkbutton(frame_Files.next,text="Next",font="our_font",width="10",
                                 command=function(){
                                   File_to_Hypo <- list(filepath_frequency, filepath_mutation, filepath_DNAtype, mode, range, as.integer(ped.max))
                                   names(File_to_Hypo) <- c("filepath_frequency", "filepath_mutation", "filepath_DNAtype", "mode", "range", "ped.max")
                                   File.result <<- checkFile(File_to_Hypo)
                                   if (class(File.result) == "list"){
                                     if (!already_pack.Linkage_Hypo || tclvalue(tkmessageBox(message="Setting of Linkage and Hypothesis will be deleted. Do you want to continue?", type="okcancel", parent = window, icon="warning")) == "ok"){
                                       if (already_pack.Linkage_Hypo){
                                         tkdestroy(frame_Linkage)
                                         tkdestroy(frame_Hypo)
                                         tkdestroy(frame_LR_result)
                                         tkdestroy(frame_LR_simu)
                                         tkdestroy(frame_Simu)
                                         already_pack.apply <<- already_pack.LR <<- already_pack.simu <<- F
                                       }
                                       else{
                                         already_pack.Linkage_Hypo <<- T
                                       }
                                       frame_Linkage <<- tkframe(tab.Linkage)
                                       frame_Hypo <<- tkframe(tab.Hypo)
                                       frame_LR_simu <<- tkframe(frame_LR)
                                       tkpack(frame_Linkage)
                                       tkpack(frame_Hypo)
                                       makeLinkage_Hypo(File.result)
                                       tk2notetab.select(Tabs, "Linkage")
                                     }
                                   }
                                                   }
                                 )

    tkpack(outframe.file,anchor="w",padx="40",pady="20")
    tkpack(frame_file.frequency,anchor="w",padx="20",pady="5")
    tkpack(frame_file.mutation,anchor="w",padx="20",pady="5")
    if (mode == "LR") tkpack(frame_file.DNAtype,anchor="w",padx="20",pady="5")

    tkpack(frame_Files.next, anchor="w", padx="40", pady="5", fill="x")

    tkpack(label_frequency,anchor="w",side="left")
    tkpack(check_frequency,anchor="w",side="left",padx="5")
    tkpack(button_frequency,padx="5")
    tkpack(label_mutation,anchor="w",side="left")
    tkpack(check_mutation,anchor="w",side="left",padx="5")
    tkpack(button_mutation,padx="5")
    if (mode == "LR"){
      tkpack(label_DNAtype,anchor="w",side="left")
      tkpack(check_DNAtype,anchor="w",side="left",padx="5")
      tkpack(button_DNAtype,padx="5")
    }
    tkpack(button_Files.next, anchor="e")
  }

  checkLinkage <- function(locus_set_name, useDefault, var_locus, linkage_number, locus1.var, locus2.var, RR.var, temp.N){
    if (is.na(as.integer(tclvalue(temp.N)))){
      tkmessageBox(message = "Please enter an integer of >2 to Number of simulation.", parent = window, icon="error")
    }
    else if (as.integer(tclvalue(temp.N)) < 2){
      tkmessageBox(message = "Please enter an integer of >2 to Number of simulation.", parent = window, icon="error")
    }
    else if (tclvalue(useDefault) == 1){
      default.sort <- sapply(default_locus_set[[File.result$locus_set_name]][[1]], function(x){return(which(x == File.result$locus_set))})
      linkage.loci <- apply(default_locus_set[[locus_set_name]][[2]], c(1,2), function(x){return(default.sort[x])})
      RR <- default_locus_set[[locus_set_name]][[3]]
      return(list(linkage.loci = linkage.loci, RR = RR, N = as.integer(tclvalue(temp.N))))
    }
    else if (tclvalue(useDefault) == -1 || linkage_number == 0){
      linkage.loci <- RR <- NA
      return(list(linkage.loci = linkage.loci, RR = RR, N = as.integer(tclvalue(temp.N))))
    }
    else{
      linkage.pair <- sapply(1:linkage_number,function(x){if(tclvalue(locus1.var[[x]]) != "" && tclvalue(locus2.var[[x]]) != "" && tclvalue(RR.var[[x]]) != "") return(x) else return(0)})
      if (is.element(0, linkage.pair)){
        tkmessageBox(message = "Error: incorrect setting of linked loci", parent = window, icon = "error")
      }
      else{
        locus1 <- sapply(linkage.pair,function(x){return(which(var_locus == tclvalue(locus1.var[[x]])))})
        locus2 <- sapply(linkage.pair,function(x){return(which(var_locus == tclvalue(locus2.var[[x]])))})
        RR <- sapply(linkage.pair,function(x){return(as.numeric(tclvalue(RR.var[[x]])))})

        if (!length(unique(c(locus1,locus2))) == length(c(locus1,locus2))){
          tkmessageBox(message = "Error: incorrect setting of linked loci", parent = window, icon="error")
        }
        else if (prod(0 <= RR) * prod(RR <= 1) == 0){
          tkmessageBox(message = "Please enter a value of 0-1 to Recombination rate.", parent = window, icon="error")
        }
        else{
          linkage.loci <- cbind(locus1,locus2)
          return(list(linkage.loci = linkage.loci, RR = RR, N = as.integer(tclvalue(temp.N))) )
        }
      }
    }
  }

  makeLinkage_Hypo <- function(File.result){
    var_locus <- File.result$locus_set
    pair.max <- floor(length(var_locus) / 2)
    
    locus1.var <- locus2.var <- RR.var <- list()
    
    combobox_Linkage.locus1 <- combobox_Linkage.locus2 <- entry_Linkage.RR <- button_Linkage.vanish <- list()

    outframe_Linkage.setting <- tkframe(frame_Linkage)
    frame_Linkage.radiobutton <- tkframe(outframe_Linkage.setting)
    frame_Linkage.add <- tkframe(outframe_Linkage.setting)
    frame_Linkage.label <- tkframe(outframe_Linkage.setting)
    frame_Linkage.select <- tkframe(outframe_Linkage.setting)
    frame_Linkage.select.input <- list()
    frame_Linkage.next <- tkframe(frame_Linkage)
    
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
        tclvalue(RR.var[[linkage_number]]) <<- "0"
        linkage_number <<- linkage_number - 1
      }
    }
    
    for (i in 1:pair.max){
      locus1.var[[i]] <- tclVar("")
      locus2.var[[i]] <- tclVar("") 
      RR.var[[i]] <- tclVar("0")
    
      frame_Linkage.select.input[[i]] <- tkframe(frame_Linkage.select)
      combobox_Linkage.locus1[[i]] <- ttkcombobox(frame_Linkage.select.input[[i]],values=var_locus,textvariable=locus1.var[[i]],font="our_font",width="10",state="disable")
      combobox_Linkage.locus2[[i]] <- ttkcombobox(frame_Linkage.select.input[[i]],values=var_locus,textvariable=locus2.var[[i]],font="our_font",width="10",state="disable")
      entry_Linkage.RR[[i]] <- tkentry(frame_Linkage.select.input[[i]], textvariable=RR.var[[i]], width="15", font="our_font", background="#ffffff" ,state="disable")
      button_Linkage.vanish[[i]] <- tkbutton(frame_Linkage.select.input[[i]],text="x",font="our_font",width=5,state="disable", command = button_vanish(i) )
    }

    cb_default.change <- function(x){
      if (x == 1){
        tkpack(frame_Linkage.add, anchor="w", padx=40, pady=5)
        tkpack(frame_Linkage.label, anchor="w", padx=40, pady=5)
        tkpack(frame_Linkage.select, padx=40, anchor="w")
        if (linkage_number > 0){
          for (i in 1:linkage_number){
            tkpack.forget(frame_Linkage.select.input[[i]])
          }
        }
        linkage_number <<- length(RR.default)
        default.sort <- sapply(default_locus_set[[File.result$locus_set_name]][[1]], function(x){return(which(x == File.result$locus_set))})
        for (i in 1:linkage_number){
          tkpack(frame_Linkage.select.input[[i]], anchor="w", pady="5")
          tkconfigure(combobox_Linkage.locus1[[i]], state="disable")
          tkconfigure(combobox_Linkage.locus2[[i]], state="disable")
          tkconfigure(entry_Linkage.RR[[i]], state="disable")
          tkconfigure(button_Linkage.vanish[[i]], state="disable")
          tclvalue(locus1.var[[i]]) <- var_locus[default.sort[linkage.loci.default[i,1]]]
          tclvalue(locus2.var[[i]]) <- var_locus[default.sort[linkage.loci.default[i,2]]]
          tclvalue(RR.var[[i]]) <- RR.default[i]
        }
        tkconfigure(button_Linkage.add, state="disable")
      }
      else if(x == 0){
        tkpack(frame_Linkage.add, anchor="w", padx=40, pady=5)
        tkpack(frame_Linkage.label, anchor="w", padx=40, pady=5)
        tkpack(frame_Linkage.select, padx=40, anchor="w")
        for (i in 1:pair.max){
          tkconfigure(combobox_Linkage.locus1[[i]], state="readonly")
          tkconfigure(combobox_Linkage.locus2[[i]], state="readonly")
          tkconfigure(entry_Linkage.RR[[i]], state="normal")
          tkconfigure(button_Linkage.vanish[[i]], state="normal")
        }
        tkconfigure(button_Linkage.add, state="normal")
      }
      else{
        tkpack.forget(frame_Linkage.add)
        tkpack.forget(frame_Linkage.label)
        tkpack.forget(frame_Linkage.select)
      }
    }
    
    button_Linkage.add <- tkbutton(frame_Linkage.add,text="Add",font="our_font",width=10,state="disable",
                                   command = function(){
                                     if (0 <= linkage_number && linkage_number < pair.max){
                                       linkage_number <<- linkage_number + 1
                                       tkpack(frame_Linkage.select.input[[linkage_number]],anchor="w", pady="5")
                                     }
                                                       }
                                   )

    if (File.result$locus_set_name != "-"){
      useDefault <- tclVar("1")
      radiobutton_off <- tkradiobutton(frame_Linkage.radiobutton, anchor="w", width=40, font="our_font", text=paste("Default settings of", File.result$locus_set_name, collapse=""), variable=useDefault, value=1, command=function(){cb_default.change(as.integer(tclvalue(useDefault)))})
      radiobutton_on <- tkradiobutton(frame_Linkage.radiobutton, anchor="w", width=40, font="our_font", text="Manual setting", variable=useDefault, value=0, command=function(){cb_default.change(as.integer(tclvalue(useDefault)))})
      radiobutton_none <- tkradiobutton(frame_Linkage.radiobutton, anchor="w", width=40, font="our_font", text="No linkage", variable=useDefault, value=-1, command=function(){cb_default.change(as.integer(tclvalue(useDefault)))})

      linkage.loci.default <- default_locus_set[[File.result$locus_set_name]][[2]]
      RR.default <- default_locus_set[[File.result$locus_set_name]][[3]]
      linkage_number <<- length(RR.default)
      default.sort <- sapply(default_locus_set[[File.result$locus_set_name]][[1]], function(x){return(which(x == File.result$locus_set))})
      for (i in 1:linkage_number){
        tkpack(frame_Linkage.select.input[[i]],anchor="w", pady=5)
        tclvalue(locus1.var[[i]]) <- var_locus[default.sort[linkage.loci.default[i,1]]]
        tclvalue(locus2.var[[i]]) <- var_locus[default.sort[linkage.loci.default[i,2]]]
        tclvalue(RR.var[[i]]) <- RR.default[i]
      }
    }
    else{
      useDefault <- tclVar("0")
      radiobutton_on <- tkradiobutton(frame_Linkage.radiobutton, anchor="w", width=40, font="our_font", text="Manual setting", variable=useDefault, value=0, command=function(){cb_default.change(as.integer(tclvalue(useDefault)))})
      radiobutton_none <- tkradiobutton(frame_Linkage.radiobutton, anchor="w", width=40, font="our_font", text="No linkage", variable=useDefault, value=-1, command=function(){cb_default.change(as.integer(tclvalue(useDefault)))})
      for (i in 1:pair.max){
        tkconfigure(combobox_Linkage.locus1[[i]], state="readonly")
        tkconfigure(combobox_Linkage.locus2[[i]], state="readonly")
        tkconfigure(entry_Linkage.RR[[i]], state="normal")
      }
      tkpack(frame_Linkage.select.input[[1]], anchor="w", pady=5)
      tkconfigure(button_Linkage.add, state="normal")
    }

    label_Linkage.locus1 <- tklabel(frame_Linkage.label,text="Locus 1",font="our_font",width="12")
    label_Linkage.locus2 <- tklabel(frame_Linkage.label,text="Locus 2",font="our_font",width="12")
    label_Linkage.RR <- tklabel(frame_Linkage.label,text="Recombination rate",font="our_font",width="15")
    
    button_Linkage.next <- tkbutton(frame_Linkage.next,text="Next",font="our_font",width="10", command = function(){tk2notetab.select(Tabs, "Hypothesis")} )
    
    tkpack(outframe_Linkage.setting, anchor="w", padx=40, pady="20")
    tkpack(frame_Linkage.radiobutton, anchor="w", padx=20, pady="15")
    tkpack(frame_Linkage.add, anchor="w", padx=40, pady=5)
    tkpack(frame_Linkage.label, anchor="w", padx=40, pady=5)
    tkpack(frame_Linkage.select, anchor="w", padx=40)
    tkpack(frame_Linkage.next, fill="both", expand="1", anchor="w", padx=10, pady=10)
    
    tkpack(label_Linkage.locus1, side="left", padx=5)
    tkpack(label_Linkage.locus2, side="left", padx=5)
    tkpack(label_Linkage.RR, anchor="w", side="left", padx=5)
    
    for (i in 1:pair.max){
      tkpack(combobox_Linkage.locus1[[i]], padx=5, side="left")
      tkpack(combobox_Linkage.locus2[[i]], padx=5, side="left")
      tkpack(entry_Linkage.RR[[i]], padx=5, side="left")
      tkpack(button_Linkage.vanish[[i]], padx=5, side="left")
    }
    
    if (File.result$locus_set_name != "-"){
      tkpack(radiobutton_off, anchor="w", padx=15)
    }
    tkpack(radiobutton_on, anchor="w", padx=15)
    tkpack(radiobutton_none, anchor="w", padx=15)

    tkpack(button_Linkage.add, padx=10, side="left")
    tkpack(button_Linkage.next, anchor="e", padx=10)

    sex.var_to_val <- function(x){
      if (tclvalue(x) == "Male"){
        return(1)
      }
      else if (tclvalue(x) == "Female"){
        return(2)
      }
      else{
        return(-1)
      }
    }
    parent.var_to_val <- function(x,var,founder, ped.number){
      ID <- which(c(var_Hypo, "") == tclvalue(var[[x]])) * (1 - as.integer(tclvalue(founder[[x]])))
      if (ID > ped.number){
        return(-1)
      }
      else{
        return(ID)
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
  
    Hypo.pedigree.check <- function(j){
      if (j == 1){
        sex <- sapply(sex.var[[1]][1:ped.number[1]],sex.var_to_val)
      }
      else{
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
        tkmessageBox(message = paste("Hypothesis", j, ": Error: incorrect setting of Sex", collapse="", sep=""), parent = window, icon="error")
      }
      else if (is.element(-1, df$pat) || is.element(-1, df$mat)){
        tkmessageBox(message = paste("Hypothesis", j, ": Error: incorrect setting of Father or Mother", collapse="", sep=""), parent = window, icon="error")
      }
      else{
        err <- try(pedigree(id = df$id, sex= df$sex, dadid = df$pat, momid = df$mat), silent = T)
        err.fc <- try(familycheck(famid=rep(1,nrow(df)), id = df$id, father.id = df$pat, mother.id = df$mat), silent = T)
        if ((class(err) == "try-error") || (class(err.fc) == "try-error")) {
          tkmessageBox(message = paste("Hypothesis", j, ": Error: incorrect setting of pedigree ", collapse="", sep=""), parent = window, icon="error")
        }
        else{
          return(df)
        }
      }
    }
  
    Hypo.pedigree <- function(df){
      if (class(df) == "data.frame"){
        window.ped <- tktoplevel(window)
        tkwm.title(window.ped, "View pedigree tree")
        unrelated <- which(df$pat == 0)[is.element(which(df$pat == 0), c(df$pat, df$mat)) == F]
        frame_ped <- tkframe(window.ped)
        tkpack(frame_ped)
        if (length(unrelated) != nrow(df)){
          labelframe_ped.related <- tk2labelframe(frame_ped, text="Pedigree tree", padding=c(20,20), labelanchor="nw", relief="groove", borderwidth="2")
          if (length(unrelated) != 0){
            pedobj <- pedigree(id = df[-unrelated,]$id, sex= df[-unrelated,]$sex, dadid = df[-unrelated,]$pat, momid = df[-unrelated,]$mat, affected = df[-unrelated,]$affected)
            pedigree_graph <- tkRplot(labelframe_ped.related,function(){plot(pedobj, id=as.character(df[-unrelated,]$names))})
          }
          else{
            pedobj <- pedigree(id = df$id, sex= df$sex, dadid = df$pat, momid = df$mat, affected = df$affected)
            pedigree_graph <- tkRplot(labelframe_ped.related,function(){ plot(pedobj, id=as.character(df$names))})
          }
          tkpack(labelframe_ped.related, padx=20, pady=10, anchor="w")
          tkpack(pedigree_graph) 
        }
        if (length(unrelated) != 0){
          labelframe_ped.unrelated <- tk2labelframe(frame_ped,text="Unrelated", labelanchor="nw", relief="groove", borderwidth="2")
          label_ped.unrelated <- tklabel(labelframe_ped.unrelated, text=paste(df[unrelated,]$names, collapse=", "), font="our_font")
          tkpack(labelframe_ped.unrelated, anchor="w", padx=20, pady=10, expand="1", fill="x")
          tkpack(label_ped.unrelated, anchor="w", padx="20") 
        }
        tkgrab.set(window.ped)
      }
    }

    if (File.result$mode == "LR"){
      var_ref <-  names(File.result[["DNAtype_matrix"]])[1:((length(File.result[["DNAtype_matrix"]]) - 1) / 2) * 2]
      var_Hypo <- c(var_ref, paste("UK",1:(ped.max - length(var_ref)),sep=""))
    }
    else{
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

    if (File.result$mode == "simu"){
      frame_Hypo.apply <- tkframe(frame_Hypo)
    }
    labelframe_Hypo <- list()
    frame_Hypo.button <- list()
    frame_Hypo.label <- list()
    for (j in 1:2){
      labelframe_Hypo[[j]] <- tk2labelframe(frame_Hypo,text=paste("Hypothesis", j, collapse="", sep=""),labelanchor="nw",relief="groove",borderwidth="2")
      frame_Hypo.button[[j]] <- tkframe(labelframe_Hypo[[j]])
      frame_Hypo.label[[j]] <- tkframe(labelframe_Hypo[[j]])
    }
    frame_Hypo.next <- tkframe(frame_Hypo)
  
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
  
    checkbutton_Hypo.change <- vector("list", length = 2)
    checkbutton_Hypo.change[[1]] <- function(){
      for (i in 1:ped.max){
        if (tclvalue(founder.var[[1]][[i]]) == 1){
          tkconfigure(combobox_Hypo.pat[[1]][[i]], state="disable")
          tkconfigure(combobox_Hypo.mat[[1]][[i]], state="disable")
          tclvalue(pat.var[[1]][[i]]) <- ""
          tclvalue(mat.var[[1]][[i]]) <- ""
        }
        else{
          tkconfigure(combobox_Hypo.pat[[1]][[i]], state="readonly")
          tkconfigure(combobox_Hypo.mat[[1]][[i]], state="readonly")
        }
      }
    }
    checkbutton_Hypo.change[[2]] <- function(){
      for (i in 1:ped.max){
        if (tclvalue(founder.var[[2]][[i]]) == 1){
          tkconfigure(combobox_Hypo.pat[[2]][[i]], state="disable")
          tkconfigure(combobox_Hypo.mat[[2]][[i]], state="disable")
          tclvalue(pat.var[[2]][[i]]) <- ""
          tclvalue(mat.var[[2]][[i]]) <- ""
        }
        else{
          tkconfigure(combobox_Hypo.pat[[2]][[i]], state="readonly")
          tkconfigure(combobox_Hypo.mat[[2]][[i]], state="readonly")
        }
      }
    }

    if (File.result$mode == "simu"){
      label_ref.num <- tklabel(frame_Hypo.apply, text="Number of known profiles", font="our_font", anchor="w")
      text_ref.num <- tkentry(frame_Hypo.apply, textvariable=ref.num, width=8, font="our_font", background="#ffffff")
      button_ref.num <- tkbutton(frame_Hypo.apply, text="Apply",font="our_font", width=8,
                                  command = function(){
                                    if (File.result$mode == "simu" && is.na(as.integer(tclvalue(ref.num)))){
                                      tkmessageBox(message = paste("Please enter an integer of 2-", File.result$ped.max, " to Number of known profiles.", sep="", collapse=""), parent = window, icon="error")
                                    }
                                    else if (File.result$mode == "simu" && (as.integer((tclvalue(ref.num))) < 2 || as.integer(tclvalue(ref.num)) > File.result$ped.max) ){
                                      tkmessageBox(message = paste("Please enter an integer of 2-", File.result$ped.max, " to Number of known profiles.", sep="", collapse=""), parent = window, icon="error")
                                    }
                                    else{
                                      if (!already_pack.apply || tclvalue(tkmessageBox(message="Setting of Hypothesis will be deleted. Do you want to continue?", type="okcancel", parent = window, icon="warning")) == "ok"){
                                        if (already_pack.apply){
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
                                        }
                                        else{
                                          already_pack.apply <<- T
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
                                    }
                                                      }
                                )
    }

    frame_Hypo.select.input <- label_Hypo.listname <- combobox_Hypo.sex <- combobox_Hypo.pat <- combobox_Hypo.mat <- checkbutton_Hypo.founder <- vector("list", length = 2)

    for (j in 1:2){
      for (i in 1:ped.max){
        frame_Hypo.select.input[[j]][[i]] <- tkframe(frame_Hypo.select[[j]])
        label_Hypo.listname[[j]][[i]] <- tklabel(frame_Hypo.select.input[[j]][[i]], text=as.list(var_Hypo)[[i]], font="our_font", width=12)
        combobox_Hypo.sex[[j]][[i]] <- ttkcombobox(frame_Hypo.select.input[[j]][[i]], values=var_sex, textvariable=sex.var[[j]][[i]], font="our_font", width=6, state="readonly")
        combobox_Hypo.pat[[j]][[i]] <- ttkcombobox(frame_Hypo.select.input[[j]][[i]], values=var_parent[[j]], textvariable=pat.var[[j]][[i]], font="our_font", width=10, state="readonly")
        combobox_Hypo.mat[[j]][[i]] <- ttkcombobox(frame_Hypo.select.input[[j]][[i]], values=var_parent[[j]], textvariable=mat.var[[j]][[i]], font="our_font", width=10, state="readonly")
        checkbutton_Hypo.founder[[j]][[i]] <- tkcheckbutton(frame_Hypo.select.input[[j]][[i]], variable=founder.var[[j]][[i]], width=5, command=checkbutton_Hypo.change[[j]])
      }
    }
    if (File.result$mode == "LR"){
      for (i in 1:length(var_ref)){
        tkconfigure(combobox_Hypo.sex[[j]][[i]], textvariable=sex.var[[1]][[i]])
      } 
    }

    button_Hypo.add <- button_Hypo.vanish <- button_Hypo.pedigree <- vector("list", length = 2)
    button_Hypo.add[[1]] <- tkbutton(frame_Hypo.button[[1]], text="Add", font="our_font", width=9, command = function(){Hypo.add(1)} )
    button_Hypo.add[[2]] <- tkbutton(frame_Hypo.button[[2]], text="Add", font="our_font", width=9, command = function(){Hypo.add(2)} )
    button_Hypo.vanish[[1]] <- tkbutton(frame_Hypo.button[[1]], text="Delete", font="our_font", width=9, command = function(){Hypo.vanish(1)} )
    button_Hypo.vanish[[2]] <- tkbutton(frame_Hypo.button[[2]], text="Delete", font="our_font", width=9, command = function(){Hypo.vanish(2)} )
    button_Hypo.pedigree[[1]] <- tkbutton(frame_Hypo.button[[1]], text="View pedigree tree",font="our_font", command = function(){Hypo.pedigree(Hypo.pedigree.check(1))} )
    button_Hypo.pedigree[[2]] <- tkbutton(frame_Hypo.button[[2]], text="View pedigree tree",font="our_font", command = function(){Hypo.pedigree(Hypo.pedigree.check(2))} )
  
    label_Hypo.name <- label_Hypo.sex <- label_Hypo.pat <- label_Hypo.mat <- label_Hypo.founder <- vector("list", length = 2)
  
    for (j in 1:2){
      label_Hypo.name[[j]] <- tklabel(frame_Hypo.label[[j]],text="Name", font="our_font", width=12)
      label_Hypo.sex[[j]] <- tklabel(frame_Hypo.label[[j]],text="Sex", font="our_font", width=8)
      label_Hypo.pat[[j]] <- tklabel(frame_Hypo.label[[j]],text="Father", font="our_font", width=12)
      label_Hypo.mat[[j]] <- tklabel(frame_Hypo.label[[j]],text="Mother", font="our_font", width=12)
      label_Hypo.founder[[j]] <- tklabel(frame_Hypo.label[[j]],text="founder",font="our_font", width=7)
    }
    if (File.result$mode == "simu"){
      label_N <- tklabel(frame_Hypo.next, text="Number of simulations", width=18, font="our_font")
      text_N <- tkentry(frame_Hypo.next, textvariable=N, width=8, font="our_font", background="#ffffff")
    }
    button_Hypo.next <- tkbutton(frame_Hypo.next, text="Calculate", font="our_font", width="10",
                                 command = function(){
                                        if (File.result$mode == "simu"){
                                          temp.N <- N
                                        }
                                        else{
                                          temp.N <- tclVar("2")
                                        }
                                        Linkage_to_calc <- checkLinkage(File.result$locus_set_name, useDefault, var_locus, linkage_number, locus1.var, locus2.var, RR.var, temp.N)
                                        if (class(Linkage_to_calc) == "list"){
                                          df <- vector("list", length = 2)
                                          df[[1]] <- Hypo.pedigree.check(1)
                                          if (class(df[[1]]) == "data.frame"){
                                            df[[2]] <- Hypo.pedigree.check(2)
                                            if(class(df[[2]]) == "data.frame"){
                                              today <- Sys.time()
                                              temp <- list(today = today, df = df)
                                              Hypo.result <<- c(File.result, Linkage_to_calc, temp)

                                              if (Hypo.result$mode == "LR"){
                                                if (!already_pack.LR || tclvalue(tkmessageBox(message="Previous results will be deleted. Do you want to continue?", type="okcancel", icon="warning")) == "ok"){
                                                  if (already_pack.LR){
                                                    tkdestroy(frame_LR_result)
                                                    tkpack.forget(frame_LR_simu)
                                                    tkdestroy(frame_Simu)
                                                    already_pack.simu <<- F
                                                  }
                                                  else{
                                                    already_pack.LR <<- T
                                                  }
                                                  LR_result <<- Calculate(Hypo.result)
                                                  frame_LR_result <<- tkframe(frame_LR)
                                                  tkpack(frame_LR_result)
                                                  tkpack(frame_LR_simu, anchor="w")
                                                  makeLR(LR_result)
                                                  tk2notetab.select(Tabs, "Result (case)")
                                                }
                                              }
                                              else{
                                                if (!already_pack.simu || tclvalue(tkmessageBox(message="Previous results will be deleted. Do you want to continue?", type="okcancel", icon="warning")) == "ok"){
                                                  if (already_pack.simu){
                                                    tkdestroy(frame_Simu)
                                                  }
                                                  else{
                                                    already_pack.simu <<- T
                                                  }
                                                  simu_result <<- Simulate(Hypo.result)
                                                  frame_Simu <<- tkframe(tab.Simu)
                                                  tkpack(frame_Simu)
                                                  makeSimu(simu_result)
                                                  tk2notetab.select(Tabs, "Result (simulation)")
                                                }
                                              }
                                            }
                                          }
                                        }
                                                     }
                                )
    if (File.result$mode == "LR"){
      tkpack(labelframe_Hypo[[1]], anchor="n", padx=20, pady=10, side="left")
      tkpack(labelframe_Hypo[[2]], anchor="n", padx=20, pady=10)
      tkpack(frame_Hypo.next, anchor="e", padx=10, pady=10)
    }
    else{
      tkpack(frame_Hypo.apply, anchor="w", padx=20, pady=20)
      tkpack(label_ref.num, anchor="w", padx=5, side="left")
      tkpack(text_ref.num, padx=5, side="left")
      tkpack(button_ref.num, padx=10)
    }

    for (j in 1:2){
      tkpack(frame_Hypo.button[[j]],anchor="w",padx="10",pady="10")
      tkpack(frame_Hypo.label[[j]],anchor="w",padx="10",pady="10")
      tkpack(frame_Hypo.select.out[[j]], anchor="w")
      if (File.result$mode == "LR"){
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
    tkconfigure(canvas_Hypo[[1]], width=505, height=192, scrollregion=c(0,0,1000,32 * length(var_ref)), yscrollcommand=function(...)tkset(scr_Hypo[[1]],...))
    tkcreate(canvas_Hypo[[1]], "window", 0, 0, anchor="nw", window=frame_Hypo.select[[1]])
    tkpack(scr_Hypo[[2]])
    tkconfigure(scr_Hypo[[2]], command=function(...)tkyview(canvas_Hypo[[2]],...))
    tkpack.configure(scr_Hypo[[2]], side="right", expand="1", fill="y")
    tkpack(canvas_Hypo[[2]], side="left")
    tkconfigure(canvas_Hypo[[2]], width=505, height=192, scrollregion=c(0,0,1000,32 * length(var_ref)), yscrollcommand=function(...)tkset(scr_Hypo[[2]],...))
    tkcreate(canvas_Hypo[[2]], "window", 0, 0, anchor="nw", window=frame_Hypo.select[[2]])

    if (File.result$mode == "simu"){
      tkpack(label_N, padx="5", pady="5", anchor="e", side="left")
      tkpack(text_N, padx="5", pady="5", anchor="e", side="left")
    }
    tkpack(button_Hypo.next, anchor="e", padx="15")

    LR.N <- tclVar("")
    labelframe_LR_simu <- tk2labelframe(frame_LR_simu, text="Simulation", labelanchor="nw", relief="groove", borderwidth=2)
    label_LR.N <- tklabel(labelframe_LR_simu, text="Number of simulations", font="our_font", width="18")
    text_LR.N <- tkentry(labelframe_LR_simu, textvariable=LR.N, width="8", font="our_font", background="#ffffff")
    button_LR.simu <- tkbutton(labelframe_LR_simu, text="Calculation", font="our_font", width="10",
                                       command=function(){
                                         Linkage_to_calc <- checkLinkage(LR_result$locus_set_name, useDefault, var_locus, linkage_number, locus1.var, locus2.var, RR.var, LR.N)
                                         if (class(Linkage_to_calc) == "list"){
                                           df <- vector("list", length = 2)
                                           df[[1]] <- Hypo.pedigree.check(1)
                                           if (class(df[[1]]) == "data.frame"){
                                             df[[2]] <- Hypo.pedigree.check(2)
                                             if(class(df[[2]]) == "data.frame"){
                                               if (!identical(list(Linkage_to_calc$linkage.loci, Linkage_to_calc$RR, df), list(LR_result$linkage.loci, LR_result$RR, LR_result$df)) ){
                                                 tkmessageBox(message="Setting of Linkage or Hypothesis has been modified.\nPlease enter Calculation button in Hypothesis tab.", icon="error")
                                               }
                                               else{
                                                 if (!already_pack.simu || tclvalue(tkmessageBox(message="Previous results will be deleted. Do you want to continue?", type="okcancel", icon="warning")) == "ok"){
                                                   if (already_pack.simu){
                                                     tkdestroy(frame_Simu)
                                                   }
                                                   else{
                                                     already_pack.simu <<- T
                                                   }
                                                   LR_result$N <- as.integer(tclvalue(LR.N))
                                                   LR_result$today <- Sys.time()
                                                   simu_result <<- Simulate(LR_result[-which(names(LR_result) == "LR")])
                                                   frame_Simu <<- tkframe(tab.Simu)
                                                   tkpack(frame_Simu)
                                                   makeSimu(simu_result)
                                                   tk2notetab.select(Tabs, "Result (simulation)")
                                                 }
                                               }
                                             }
                                           }
                                         }
                                                         }
                               )
    tkpack(labelframe_LR_simu, anchor="w", padx=40, pady=10)
    tkpack(label_LR.N, anchor="w", side="left")
    tkpack(text_LR.N, anchor="w", side="left", padx=10)
    tkpack(button_LR.simu, anchor="w", padx=10, pady=5)
  }

  makeLR <- function(LR_result){
    if (!is.na(LR_result$RR[1])){
      linkage.loci.vector <- as.vector(LR_result$linkage.loci)
      log_allelename <- paste(paste(LR_result$locus_set[-linkage.loci.vector], ": "), collapse="\n")
      log_LR <- paste(signif(LR_result$LR[[1]],digit=3), collapse="\n")
      LRprod <- signif(prod(LR_result$LR[[1]]) * prod(LR_result$LR[[2]]),digit=3)
      log_LR.locus <- tclVar(paste("Overall LR : \n\n", log_allelename, sep="", collapse=""))
      log_LR.value <- tclVar(paste(LRprod, "\n\n", log_LR, sep="", collapse=""))

      log_linkage.loci <- apply(LR_result$linkage.loci,c(1,2),function(x){return(LR_result$locus_set[x])})
      log_LR.linkage.locus <- tclVar(paste("\n\n", paste(paste(log_linkage.loci[,1], " - ", log_linkage.loci[,2], ": "), collapse="\n"),sep="",collapse=""))
      log_LR.linkage.value <- tclVar(paste("\n\n", paste(signif(LR_result$LR[[2]],digit=3), collapse="\n"),sep="",collapse=""))
    }
    else{
      log_allelename <- paste(paste(LR_result$locus_set, ": "), collapse="\n")
      log_LR <- paste(signif(LR_result$LR,digit=3), collapse="\n")
      LRprod <- signif(prod(LR_result$LR),digit=3)
      log_LR.locus <- tclVar(paste("Overall LR : \n\n", log_allelename, sep="", collapse=""))
      log_LR.value  <- tclVar(paste(LRprod, "\n\n", log_LR, sep="", collapse=""))
      log_LR.linkage.locus <- tclVar("\n\nNo linkage")
      log_LR.linkage.value <- tclVar("")
    }

    labelframe_LR <- tk2labelframe(frame_LR_result,text="Likelihood ratio",labelanchor="nw",relief="groove",borderwidth="2")
    frame_LR.locus.value <- tkframe(labelframe_LR)
    frame_LR.linkage <- tkframe(labelframe_LR)
    frame_LR.linkage.locus.value <- tkframe(frame_LR.linkage)
    frame_LR.report <- tkframe(frame_LR.linkage)
  
    label_LR.locus <- tklabel(frame_LR.locus.value, textvariable=log_LR.locus, font="our_font", justify="left")
    label_LR.value <- tklabel(frame_LR.locus.value, textvariable=log_LR.value, font="our_font", justify="left")
    label_LR.linkage.locus <- tklabel(frame_LR.linkage.locus.value, textvariable=log_LR.linkage.locus, font="our_font", justify="left")
    label_LR.linkage.value <- tklabel(frame_LR.linkage.locus.value, textvariable=log_LR.linkage.value, font="our_font", justify="left")
    button_LR.report <- tkbutton(frame_LR.report, text="Report", font="our_font", command=function() LR_Report(LR_result))
  
    tkpack(labelframe_LR, anchor="w", padx=40, pady=10)
    tkpack(frame_LR.locus.value, anchor="w", side="left")
    tkpack(frame_LR.linkage, fill="both", anchor="w", padx="20", expand="1", pady="15")
    tkpack(frame_LR.linkage.locus.value, fill="none", anchor="w")
    tkpack(frame_LR.report,fill="x", anchor="s", expand="1", pady="10")

    tkpack(label_LR.locus, anchor="w",side="left", padx="10", pady="15")
    tkpack(label_LR.value, anchor="w",side="left", padx="10", pady="15")
    tkpack(label_LR.linkage.locus, anchor="w", side="left")
    tkpack(label_LR.linkage.value, anchor="w", padx="10")
    tkpack(button_LR.report,fill="x")
  }

  makeSimu <- function(simu_result){
    simu_items <- c("","min","5th percentile","25th percentile","50th percentile","75th percentile","95th percentile","max")
    simu_title <- list()

    simu_value <- vector("list", length = 2)
    LRprod <- list()
    for (k in 1:2){
      if (!is.na(simu_result$RR[1])){
        LRprod[[k]] <- sort(apply(simu_result$LR[[k]][[1]],1,prod) * apply(simu_result$LR[[k]][[2]],1,prod))
      }
      else{
        LRprod[[k]] <- sort(apply(simu_result$LR[[k]],1,prod))
      }
    }

    log_simu <- list(c(), c())
    for (k in 1:2){
      log_simu[[k]][1] <- paste("H", k, " true testing", sep="", collapse="")
      log_simu[[k]][2] <- signif(quantile(LRprod[[k]])[["0%"]], digit=3)
      log_simu[[k]][3] <- signif(LRprod[[k]][floor(length(LRprod[[k]]) / 20) + 1], digit=3)
      log_simu[[k]][4] <- signif(quantile(LRprod[[k]])["25%"], digit=3)
      log_simu[[k]][5] <- signif(quantile(LRprod[[k]])["50%"], digit=3)
      log_simu[[k]][6] <- signif(quantile(LRprod[[k]])["75%"], digit=3)
      log_simu[[k]][7] <- signif(LRprod[[k]][ceiling(length(LRprod[[k]]) * 19 / 20)], digit=3)
      log_simu[[k]][8] <- signif(quantile(LRprod[[k]])["100%"], digit=3)
    }
 
    outframe_simu <- tkframe(frame_Simu)
    frame_simu.result <- tkframe(outframe_simu)
    frame_simu.button <- tkframe(outframe_simu)
  
    label_simu.title <- tklabel(frame_simu.result, width=15, justify="right", text=paste(simu_items, sep="", collapse="\n"), font="our_font")
    label_simu.value <- list()
    for (k in 1:2){
      label_simu.value[[k]] <- tklabel(frame_simu.result, width=15, text=paste(log_simu[[k]], sep="", collapse="\n"), font="our_font")
    }
  
    button_simu.graph <- tkbutton(frame_simu.button,text=" Graph ", font="our_font",
                                          command=function(){
                                            xmax <- max(c(density(log10(LRprod[[1]]))[["x"]], density(log10(LRprod[[2]]))[["x"]]))
                                            xmin <- min(c(density(log10(LRprod[[1]]))[["x"]], density(log10(LRprod[[2]]))[["x"]]))
                                            xwidth <- xmax - xmin
                                            xmax <- xmax + xwidth / 15
                                            xmin <- xmin - xwidth / 15
                                            ymax <- max(c(density(log10(LRprod[[1]]))[["y"]], density(log10(LRprod[[2]]))[["y"]])) * 1.1
                                            graph.draw <- function(){
                                                            plot(density(log10(LRprod[[1]])), col=rgb(1,0,0), lwd=2, cex.lab=1.5, xlim=c(xmin,xmax), ylim=c(0,ymax), xlab="", ylab="", ann=F, main="") 
                                                            par(new=T)
                                                            plot(density(log10(LRprod[[2]])), col=rgb(0,0.5,1), lwd=2, cex.lab=1.5, xlim=c(xmin,xmax), ylim=c(0,ymax), xlab=expression(log[10](LR)), main="")
                                                            legend("topright",col=c(rgb(1,0,0),rgb(0,0.5,1)),lty=c(1,1),lwd=c(3,3),c("H1 true","H2 true"))
                                                          }
                                            window.graph <- tktoplevel(window)
                                            tkwm.title(window.graph, "Graph")
                                            frame.graph <- tkframe(window.graph, padx=20, bg="#ffffff")
                                            simu_graph <- tkRplot(frame.graph, graph.draw)
                                            tkpack(frame.graph)
                                            tkpack(simu_graph)
                                            tkgrab.set(window.graph)
                                          }
                                        )
    button_simu.report <- tkbutton(frame_simu.button,text=" Report ", font="our_font", command=function() simu_Report(simu_result))
    button_simu.download <- tkbutton(frame_simu.button,text=" All data ", font="our_font", command=function() alldata_Report(simu_result))
  
    tkpack(outframe_simu, anchor="w", padx=20, pady=40)
    tkpack(frame_simu.result, anchor="w")
    tkpack(frame_simu.button, anchor="w", padx=10, pady=20)

    tkpack(label_simu.title, anchor="w", side="left")
    for (k in 1:2){
      tkpack(label_simu.value[[k]], anchor="w", side="left")
    }
    tkpack(button_simu.graph, anchor="w", side="left", padx=10)
    tkpack(button_simu.report, anchor="w", side="left", padx=10)
    tkpack(button_simu.download, anchor="w", padx=10)
  }

  pedigree.simple <- function(df){
    df.simple <- pedigree.simple.loop(df)
    if (is.element("x", df.simple$affected)){
      vanish_number <- which(df.simple$affected == "x")
      df.simple <- df.simple[-vanish_number, ]
    }
    return(df.simple)
  }

  pedigree.simple.loop <- function(df){
    banish_number <- c()
    for (i in 1:nrow(df)){
      if ((df$affected[i] != "x") && (df$DNAtype[i] == 0) && (!is.element(df$id[i], c(df$pat, df$mat))) ){
        banish_number <- c(banish_number,i)
      } 
      else if ((length(which(df$pat[i] == c(df$pat, df$mat))) <= 1) &&
               (length(which(df$mat[i] == c(df$pat, df$mat))) <= 1) &&
                df[df$id == df$pat[i],]$DNAtype == 0 && 
                df[df$id == df$mat[i],]$DNAtype == 0 &&
                df[df$id == df$pat[i],]$pat == 0 &&
                df[df$id == df$mat[i],]$mat == 0 ){
        df$pat[i] <- 0
        df$mat[i] <- 0
        df <- pedigree.simple.loop(df)
      }
    }
    if (!is.null(banish_number)){
      b <- min(banish_number)
      df$affected[b] <- "x"
      df$pat[b] <- df$mat[b] <- 0
      df <- pedigree.simple.loop(df)
    }
    return(df)
  }

  check.pedigree <- function(df1,df2){
    res <- 1
    if(nrow(df1) != nrow(df2)){
      res <- 0
    }else if(length(which(df1[,3]==1))!=length(which(df2[,3]==1))){
      res <- 0
    }else{
      Ref <- df1$DNAtype[which(df1$DNAtype!=0)]
      Ref.list <- list()
      Ref.list[[1]] <- list()
      Ref.list[[2]] <- list()
      for(i in Ref){
        ance.H1 <- list()
        ance.H1[[1]] <- matrix(as.numeric(df1[which(df1$DNAtype==i),c(4,5)]),1,2)
        ance.H2 <- list()
        ance.H2[[1]] <- matrix(as.numeric(df2[which(df2$DNAtype==i),c(4,5)]),1,2)
        for(j in 2:5){
          ance.H1[[j]] <- matrix(0,2*nrow(ance.H1[[j-1]]),2)
          ance.H2[[j]] <- matrix(0,2*nrow(ance.H1[[j-1]]),2)
          for(k in 1:nrow(ance.H1[[j-1]])){
            fa.id.H1 <- ance.H1[[j-1]][k,1]
            if(fa.id.H1!=0){
              ance.H1[[j]][2*k-1,] <- as.numeric(df1[which(df1$id==fa.id.H1),c(4,5)])
            }
            mo.id.H1 <- ance.H1[[j-1]][k,2]
            if(mo.id.H1!=0){
              ance.H1[[j]][2*k,] <- as.numeric(df1[which(df1$id==mo.id.H1),c(4,5)])
            }
            fa.id.H2 <- ance.H2[[j-1]][k,1]
            if(fa.id.H2!=0){
              ance.H2[[j]][2*k-1,] <- as.numeric(df2[which(df2$id==fa.id.H2),c(4,5)])
            }
            mo.id.H2 <- ance.H2[[j-1]][k,2]
            if(mo.id.H2!=0){
              ance.H2[[j]][2*k,] <- as.numeric(df2[which(df2$id==mo.id.H2),c(4,5)])
            }
          }
          if(all(c(as.vector(ance.H1[[j]]),as.vector(ance.H2[[j]]))==0)){
            break
          }
        }
        Ref.list[[1]][[length(Ref.list[[1]])+1]] <- unlist(ance.H1)
        Ref.list[[2]][[length(Ref.list[[2]])+1]] <- unlist(ance.H2)      
        if(!identical(which(Ref.list[[1]][[length(Ref.list[[1]])]]==0),which(Ref.list[[2]][[length(Ref.list[[2]])]]==0))){
          res <- 0
        }
      }
   
      Ref.H1 <- unlist(Ref.list[[1]])
      Ref.H2 <- unlist(Ref.list[[2]])
      Ref.H1.nof <- Ref.H1[-which(Ref.H1==0)]
      Ref.H2.nof <- Ref.H1[-which(Ref.H2==0)]
      for(i in 1:length(Ref.H1.nof)){
        id.H1 <- Ref.H1.nof[i]
        id.pos.H1 <- which(Ref.H1.nof==id.H1)
        id.H2 <- Ref.H2.nof[i]
        id.pos.H2 <- which(Ref.H2.nof==id.H2)
        if(!identical(id.pos.H1,id.pos.H2)){ 
          res <- 0
          break
        }
      }
      }
    return(res)
  }

  pedigree.numbering <- function(x,names){
    if(x == "0"){
      return(0)
    }
    else{
      return(which(names == x))
    }
  }

  PES_function <- function(calc.info){
    pedigree <- apply(calc.info$df[,c("pat","mat")], c(1,2), pedigree.numbering, calc.info$df[,"id"])
    DNAtype <- calc.info$df[,"DNAtype"]
    alleletype <- calc.info$alleletype
    target_number <- calc.info$target_number
    pa_mutation <- calc.info$mutation_pat
    ma_mutation <- calc.info$mutation_mat
    al <- calc.info$aList
    pl <- calc.info$pList
    range <- 0
  
    result.LR <- matrix(1,dim(alleletype)[4],length(al))
   
    n <- dim(pedigree)[1]
    founder <- which(pedigree[,1] == 0) 
    non_founder <- which(pedigree[,1] != 0)
    known_numbers <- which(calc.info$df$DNAtype != 0)
  
    if (length(non_founder) >= 1){ 
      nodekind <- c("pg","mg","gt")
      node.person <- unlist(lapply(paste("n",1:n,sep=""),function(x){return(paste(x,nodekind,sep=""))}))
      node.target <- c(paste("p",target_number,c("pg","mg"),sep=""),"target")
      node.mutation <- unlist(lapply(paste("o",non_founder,sep=""),function(x){return(paste(x,c("pg","mg"),sep=""))}))
      nodenames <- c(node.person,node.target,node.mutation)
  
      PES <- empty.graph(nodenames)
  
      arc.putative_or_true <- rbind(c(node.person[3 * target_number - 2], node.target[1]),
                                    c("target", node.target[1]),
                                    c(node.person[3 * target_number - 1], node.target[2]),
                                    c("target", node.target[2]))
  
      node.person2 <- node.person
      node.person2[(3 * target_number - 2):(3 * target_number - 1)] <- node.target[1:2]
      arc.allele_to_genotype <- matrix(nrow=2 * n,ncol=2)
      arc.allele_to_genotype[,1] <- as.vector(sapply(1:n,function(x){return(c(node.person2[3 * x - 2],node.person2[3 * x - 1]))}))
      arc.allele_to_genotype[,2] <- as.vector(sapply(1:(2 * n),function(x){return(node.person2[3 * ceiling(x / 2)])}))
  
      arc.parent_to_child <- matrix(nrow=4 * length(non_founder),ncol=2)
      for (i in 1:length(non_founder)){
        arc.parent_to_child[(4 * i - 3):(4 * i),] <- rbind(c(node.person[3 * pedigree[non_founder[i],1] - 2],node.mutation[2 * i - 1]),
                                                           c(node.person[3 * pedigree[non_founder[i],1] - 1],node.mutation[2 * i - 1]),
                                                           c(node.person[3 * pedigree[non_founder[i],2] - 2],node.mutation[2 * i]),
                                                           c(node.person[3 * pedigree[non_founder[i],2] - 1],node.mutation[2 * i]))
      }
  
      arc.mutation <- matrix(nrow=2 * length(non_founder),ncol=2)
      arc.mutation[,1] <- node.mutation
      arc.mutation[,2] <- as.vector(sapply(non_founder,function(x){return(c(node.person[3 * x - 2],node.person[3 * x - 1]))}))
  
      arc.set <- rbind(arc.putative_or_true, arc.allele_to_genotype, arc.parent_to_child, arc.mutation)
      arcs(PES) <- arc.set
  
      for(m in 1:(dim(alleletype)[4])){  
        for (k in 1:length(al)){
    
          known_gt <- apply(alleletype[,k,DNAtype,m,drop=FALSE],3,function(x){return(paste(min(x),"-",max(x),sep=""))})
          aList_number <- sort(unique(as.vector(alleletype[,k,,m])))
          uni_aList <- unique(c(al[[k]],aList_number))
          range_aList <- sort(unique(as.vector(sapply(c(-range:range),function(x){return(x + uni_aList)}))))
          aList <- range_aList[which(min(uni_aList) <= range_aList & range_aList <= max(uni_aList))]
        
          cp_convert <- function(x){
            if (any(al[[k]] == x)){
              return((pl[[k]][which(al[[k]] == x)] + 1) / (sum(pl[[k]]) + length(aList)))
            }
            else{
              return(1 / (sum(pl[[k]]) + length(aList)))
            }
          }    
          cpList <- sapply(aList,cp_convert)
    
          expected_allele <- sort(unique(as.vector(sapply(c(-range:range),function(x){return(x + aList_number)}))))
          aList2 <- c(expected_allele,99)
          cpList2 <- c(cpList[which(aList %in% expected_allele)],sum(cpList[-which(aList %in% expected_allele)]))
          pol <- length(aList2)
          gList <- unlist(lapply(aList2,function(x){return(paste(x,"-",aList2,sep=""))}))
          del <- which((1:pol^2 - 1) %/% pol > (1:pol^2 - 1) %% pol)
          gList2 <- gList[-del]
          Boolean <- c("T", "F")
    
          founder_matrix <- matrix(cpList2, ncol = length(aList2), dimnames = list(NULL, aList2))
    
          target_query <- matrix(c(0.5,0.5), ncol = 2, dimnames = list(NULL, Boolean))
    
          putative_or_true <- array(0,dim=c(pol,pol,2))
          putative_or_true[,,1] <- diag(pol)
          putative_or_true[,,2] <- cpList2
          dimnames(putative_or_true) <- list(aList2,aList2,Boolean)
    
          allele_to_genotype <- array(0,dim=c(pol ^ 2,pol,pol))
          for (i in 1:pol){
            allele_to_genotype[((i - 1) * pol + 1):(i * pol),,i] <- diag(pol)
            allele_to_genotype[(1:pol - 1) * pol + i,,i] <- diag(pol)
          }
          allele_to_genotype <- allele_to_genotype[-del,,]
          dimnames(allele_to_genotype) <- list(gList2,aList2,aList2)
    
          parent_to_child <- array(1:pol^3 - 1,dim=c(pol,pol,pol))
          parent_to_child <- apply(parent_to_child,c(1,2,3),function(x){return(sum(c((x %/% pol) %% pol + 1,x %/% (pol ^ 2) + 1) == x %% pol + 1) / 2)})
          dimnames(parent_to_child) <- list(aList2,aList2,aList2)
    
          pre_mutation_paternal <- matrix(0,nrow=pol,ncol=pol)
            for (i in 1:length(aList2)){
              for (g in -2:2){
                if (any(aList2 == aList2[i] + g)){
                  pre_mutation_paternal[i,which(aList2 == aList2[i] + g)] <- pa_mutation[k,(g+3)]
                }
              }
            }
          mutation_paternal <- t(apply(pre_mutation_paternal,1,function(x){return(x / sum(x))}))
          dimnames(mutation_paternal) <- list(aList2,aList2)
        
          pre_mutation_maternal <- matrix(0,nrow=pol,ncol=pol)
          for (i in 1:length(aList2)){
            for (g in -2:2){
              if (any(aList2 == aList2[i] + g)){
                pre_mutation_maternal[i,which(aList2 == aList2[i] + g)] <- ma_mutation[k,(g+3)]
              }
            }
          }
          mutation_maternal <- t(apply(pre_mutation_maternal,1,function(x){return(x / sum(x))}))
          dimnames(mutation_maternal) <- list(aList2,aList2)
          matrix.person <- list()
          matrix.target <- list()
          matrix.mutation <- list()
    
          for (i in founder){
            matrix.person[[3 * i - 2]] <- founder_matrix
            matrix.person[[3 * i - 1]] <- founder_matrix 
          }
    
          matrix.target[[3]] <- target_query
    
          matrix.target[[1]] <- putative_or_true
          names(dimnames(matrix.target[[1]])) <- list(node.target[1],node.person[3 * target_number - 2],"target")
          matrix.target[[2]] <- putative_or_true
          names(dimnames(matrix.target[[2]])) <- list(node.target[2],node.person[3 * target_number - 1],"target")
    
          for (i in 1:n){
            matrix.person[[3 * i]] <- allele_to_genotype
            names(dimnames(matrix.person[[3 * i]])) <- list(node.person2[3 * i],node.person2[3 * i - 2],node.person2[3 * i - 1])
          }
    
          for (i in 1:length(non_founder)){
            matrix.mutation[[2 * i - 1]] <- parent_to_child
            names(dimnames(matrix.mutation[[2 * i - 1]])) <- list(node.mutation[2 * i - 1],node.person[3 * pedigree[non_founder[i],1] - 2],node.person[3 * pedigree[non_founder[i],1] - 1])
            matrix.mutation[[2 * i]] <- parent_to_child
            names(dimnames(matrix.mutation[[2 * i]])) <- list(node.mutation[2 * i],node.person[3 * pedigree[non_founder[i],2] - 2],node.person[3 * pedigree[non_founder[i],2] - 1])
          }
    
          for (i in 1:length(non_founder)){
            matrix.person[[3 * non_founder[i] - 2]] <- mutation_paternal
            names(dimnames(matrix.person[[3 * non_founder[i] - 2]])) <- list(node.person[3 * non_founder[i] - 2],node.mutation[2 * i - 1])
            matrix.person[[3 * non_founder[i] - 1]] <- mutation_maternal
            names(dimnames(matrix.person[[3 * non_founder[i] - 1]])) <- list(node.person[3 * non_founder[i] - 1],node.mutation[2 * i])
          }
          matrix_list <- c(matrix.person,matrix.target,matrix.mutation)
    
          names(matrix_list) <- nodenames
          dfit <- custom.fit(PES,dist = matrix_list) 
          fitted.grain <-  as.grain(dfit)
          known_gtnodes <- c(nodenames[3 * known_numbers])
          net2 <- setEvidence(fitted.grain, nodes=known_gtnodes, states=known_gt) 
          result.part <- querygrain(net2,nodes=c("target"))
          LR.locus <- result.part$target["T"] / result.part$target["F"]
          result.LR[m,k] <- LR.locus
        }
      }
    }
    return(result.LR)
  }

  PES_function_linkage <- function(calc.info){
    pedigree <- apply(calc.info$df[,c("pat","mat")], c(1,2), pedigree.numbering, calc.info$df[,"id"])
    DNAtype <- calc.info$df[,"DNAtype"]
    alleletype <- calc.info$alleletype
    target_number <- calc.info$target_number
    pa_mutation <- calc.info$mutation_pat
    ma_mutation <- calc.info$mutation_mat
    al <- calc.info$aList
    pl <- calc.info$pList
    range <- 0
    linkage.loci <- calc.info$linkage.loci
    RR <- calc.info$RR
  
    result.LR <- list()
    result.LR.noL <- matrix(NA,dim(alleletype)[4],length(al))
    result.LR.L <- matrix(1,dim(alleletype)[4],nrow(linkage.loci))
  
    linkage.loci.name <- sapply(1:nrow(linkage.loci), function(i){paste(linkage.loci[i,1],"-",linkage.loci[i,2],sep="")})
    colnames(result.LR.L) <- linkage.loci.name
  
    calc.info.noL <- calc.info 
    linkage.loci.vector <- as.vector(calc.info$linkage.loci)
    calc.info.noL$aList <- calc.info$aList[-linkage.loci.vector]
    calc.info.noL$pList <- calc.info$pList[-linkage.loci.vector]
    calc.info.noL$mutation_pat <- calc.info$mutation_pat[-linkage.loci.vector, ]
    calc.info.noL$mutation_mat <- calc.info$mutation_mat[-linkage.loci.vector, ]
    calc.info.noL$alleletype <- calc.info$alleletype[,-linkage.loci.vector, , , drop=F]
    result.LR.noL <- PES_function(calc.info.noL)
  
    n <- dim(pedigree)[1]
    founder <- which(pedigree[,1] == 0)
    non_founder <- which(pedigree[,1] != 0)
    known_numbers <- which(calc.info$df$DNAtype != 0)
  
    if (length(non_founder) >= 1){ 
      nodekind <- c("pg","mg","gt")
      node.person <- unlist(lapply(paste("n",1:n,sep=""),function(x){return(paste(x,nodekind,sep=""))}))
      node.target <- c(paste("p",target_number,c("pg","mg"),sep=""))
      node.mutation <- unlist(lapply(paste("o",non_founder,sep=""),function(x){return(paste(x,c("pg","mg"),sep=""))}))
      node.linkage <- unlist(lapply(paste("o",non_founder,sep=""),function(x){return(paste(x,c("pg_fp","mg_fp"),sep=""))}))
  
      nodenames <- c(node.person,node.target,node.mutation,node.linkage)
      nodenames_L1 <- c(paste(nodenames,"_L1",sep=""))
      nodenames_L2 <- c(paste(nodenames,"_L2",sep=""))
  
      PES <- empty.graph(c(nodenames_L1,nodenames_L2,"target"))
    
      node.person2 <- node.person
      node.person2[(3 * target_number - 2):(3 * target_number - 1)] <- node.target[1:2]
      arc.allele_to_genotype <- matrix(nrow=2 * n,ncol=2)
      arc.allele_to_genotype[,1] <- as.vector(sapply(1:n,function(x){return(c(node.person2[3 * x - 2],node.person2[3 * x - 1]))}))
      arc.allele_to_genotype[,2] <- as.vector(sapply(1:(2 * n),function(x){return(node.person2[3 * ceiling(x / 2)])}))
    
      arc.parent_to_child <- matrix(nrow=6 * length(non_founder),ncol=2)
      for (i in 1:length(non_founder)){
        arc.parent_to_child[(6 * i - 5),] <- c(node.person[3 * pedigree[non_founder[i],1] - 2],node.mutation[2 * i - 1])
        arc.parent_to_child[(6 * i - 4),] <- c(node.person[3 * pedigree[non_founder[i],1] - 1],node.mutation[2 * i - 1])
        arc.parent_to_child[(6 * i - 3),] <- c(node.linkage[2 * i - 1],node.mutation[2 * i - 1])
        arc.parent_to_child[(6 * i - 2),] <- c(node.person[3 * pedigree[non_founder[i],2] - 2],node.mutation[2 * i])
        arc.parent_to_child[(6 * i - 1),] <- c(node.person[3 * pedigree[non_founder[i],2] - 1],node.mutation[2 * i])
        arc.parent_to_child[(6 * i),] <- c(node.linkage[2 * i],node.mutation[2 * i])
      }
    
      arc.L1_to_L2 <- matrix(nrow=2 * length(non_founder),ncol=2)
      arc.L1_to_L2[,1] <- c(paste(node.linkage,"_L1",sep=""))
      arc.L1_to_L2[,2] <- c(paste(node.linkage,"_L2",sep=""))
    
      arc.mutation <- matrix(nrow=2 * length(non_founder),ncol=2)
      arc.mutation[,1] <- node.mutation
      arc.mutation[,2] <- as.vector(sapply(non_founder,function(x){return(c(node.person[3 * x - 2],node.person[3 * x - 1]))}))
    
      arc.set <- matrix(nrow=8*length(non_founder)+2*n,ncol=2)
      arc.set[c(1:(2*n)),] <- arc.allele_to_genotype
      arc.set[c((2*n+1):(2*n+6*length(non_founder))),] <- arc.parent_to_child
      arc.set[c((6*length(non_founder)+2*n+1):(8*length(non_founder)+2*n)),] <- arc.mutation
    
      arc.set_L1 <- matrix(c(paste(arc.set,"_L1",sep="")),nrow(arc.set),2)
      arc.set_L2 <- matrix(c(paste(arc.set,"_L2",sep="")),nrow(arc.set),2)
    
      arc.putative_or_true <- matrix(nrow=4*ncol(linkage.loci),ncol=2)
      arc.putative_or_true[1,] <- c(paste(node.person[3 * target_number - 2],"_L1",sep=""),paste(node.target[1],"_L1",sep=""))
      arc.putative_or_true[2,] <- c("target", paste(node.target[1],"_L1",sep=""))
      arc.putative_or_true[3,] <- c(paste(node.person[3 * target_number - 1],"_L1",sep=""), paste(node.target[2],"_L1",sep=""))
      arc.putative_or_true[4,] <- c("target", paste(node.target[2],"_L1",sep=""))
      arc.putative_or_true[5,] <- c(paste(node.person[3 * target_number - 2],"_L2",sep=""),paste(node.target[1],"_L2",sep=""))
      arc.putative_or_true[6,] <- c("target", paste(node.target[1],"_L2",sep=""))
      arc.putative_or_true[7,] <- c(paste(node.person[3 * target_number - 1],"_L2",sep=""), paste(node.target[2],"_L2",sep=""))
      arc.putative_or_true[8,] <- c("target", paste(node.target[2],"_L2",sep=""))
      
      arcs_PES <- matrix(nrow=(4*ncol(linkage.loci)+18*length(non_founder)+4*n),ncol=2)
      arcs_PES[c(1:(4*ncol(linkage.loci))),] <- arc.putative_or_true
      arcs_PES[c((4*ncol(linkage.loci)+1):(4*ncol(linkage.loci)+nrow(arc.set))),] <- arc.set_L1
      arcs_PES[c((4*ncol(linkage.loci)+nrow(arc.set)+1):(4*ncol(linkage.loci)+2*nrow(arc.set))),] <- arc.set_L2
      arcs_PES[c((4*ncol(linkage.loci)+2*nrow(arc.set)+1):(4*ncol(linkage.loci)+18*length(non_founder)+4*n)),] <- arc.L1_to_L2
      arcs(PES) <- arcs_PES 
    
      matrix.list_L1 <- list()
      matrix.list_L2 <- list()
        
      for(m in 1:dim(alleletype)[4]){
        for (l in 1:nrow(linkage.loci)){
          for (k in linkage.loci[l,]){
            aList_number <- sort(unique(as.vector(alleletype[,k,,m])))
            uni_aList <- unique(c(al[[k]],aList_number))
            aList <- sort(unique(as.vector(sapply(c(-range:range),function(x){return(x + uni_aList)}))))
        
            cp_convert <- function(x){
              if (any(al[[k]] == x)){
                return((pl[[k]][which(al[[k]] == x)] + 1) / (sum(pl[[k]]) + length(aList)))
              }else{
                return(1 / (sum(pl[[k]]) + length(aList)))
              }
            }    
            cpList <- sapply(aList,cp_convert)
    
            expected_allele <- sort(unique(as.vector(sapply(c(-range:range),function(x){return(x + aList_number)}))))
            aList2 <- c(expected_allele,99)
            cpList2 <- c(cpList[which(aList %in% expected_allele)],sum(cpList[-which(aList %in% expected_allele)]))
            pol <- length(aList2)
            gList <- unlist(lapply(aList2,function(x){return(paste(x,"-",aList2,sep=""))}))
            del <- which((1:pol^2 - 1) %/% pol > (1:pol^2 - 1) %% pol)
            gList2 <- gList[-del]
            Boolean <- c("T", "F")
    
            founder_matrix <- matrix(cpList2, ncol = length(aList2), dimnames = list(NULL, aList2))
    
            target_query <- matrix(c(0.5,0.5), ncol = 2, dimnames = list(NULL, Boolean))
        
            L1_p_or_m <- matrix(c(0.5,0.5), ncol = 2, dimnames = list(NULL, Boolean))
    
            putative_or_true <- array(0,dim=c(pol,pol,2))
            putative_or_true[,,1] <- diag(pol)
            putative_or_true[,,2] <- cpList2
            dimnames(putative_or_true) <- list(aList2,aList2,Boolean)
        
            allele_to_genotype <- array(0,dim=c(pol ^ 2,pol,pol))
            for (i in 1:pol){
              allele_to_genotype[((i - 1) * pol + 1):(i * pol),,i] <- diag(pol)
              allele_to_genotype[(1:pol - 1) * pol + i,,i] <- diag(pol)
            }
            allele_to_genotype <- allele_to_genotype[-del,,]
            dimnames(allele_to_genotype) <- list(gList2,aList2,aList2)
    
            parent_to_child <- array(0,dim=c(pol,pol,pol,2))
            for(i in 1:pol){
              parent_to_child[i,,i,1] <- rep(1,pol)
              parent_to_child[i,i,,2] <- rep(1,pol)
            }
            dimnames(parent_to_child) <- list(aList2,aList2,aList2,Boolean)
    
            pre_mutation_paternal <- matrix(0,nrow=pol,ncol=pol)
              for (i in 1:length(aList2)){
                for (g in -2:2){
                  if (any(aList2 == aList2[i] + g)){
                    pre_mutation_paternal[i,which(aList2 == aList2[i] + g)] <- pa_mutation[k,(g+3)]
                  }
                }
              }
            mutation_paternal <- t(apply(pre_mutation_paternal,1,function(x){return(x / sum(x))}))
            dimnames(mutation_paternal) <- list(aList2,aList2)
        
            pre_mutation_maternal <- matrix(0,nrow=pol,ncol=pol)
            for (i in 1:length(aList2)){
              for (g in -2:2){
                if (any(aList2 == aList2[i] + g)){
                  pre_mutation_maternal[i,which(aList2 == aList2[i] + g)] <- ma_mutation[k,(g+3)]
                }
              }
            }
            mutation_maternal <- t(apply(pre_mutation_maternal,1,function(x){return(x / sum(x))}))
            dimnames(mutation_maternal) <- list(aList2,aList2)
            matrix.person <- list()
            matrix.target <- list()   
            matrix.parent_to_child <- list()
            matrix.linkage <- list()
    
            L1_L2 <- matrix(RR[l],nrow=2,ncol=2)
            diag(L1_L2) <- 1-RR[l]
            dimnames(L1_L2) <- list(Boolean,Boolean)
    
            matrix.target[[1]] <- putative_or_true
            matrix.target[[2]] <- putative_or_true
       
            if(k==linkage.loci[l,1]){
            
              names(dimnames(matrix.target[[1]])) <- list(paste(node.target[1],"_L1",sep=""),paste(node.person[3 * target_number - 2],"_L1",sep=""),"target")
              names(dimnames(matrix.target[[2]])) <- list(paste(node.target[2],"_L1",sep=""),paste(node.person[3 * target_number - 1],"_L1",sep=""),"target")
          
              for (i in founder){
                matrix.person[[3 * i - 2]] <- matrix(cpList2, ncol = length(aList2), dimnames = list(paste(node.person2[3*i-2],"_L1",sep=""), aList2))
                matrix.person[[3 * i - 1]] <- matrix(cpList2, ncol = length(aList2), dimnames = list(paste(node.person2[3*i-1],"_L1",sep=""), aList2))
              }
              for (i in non_founder){
                I <- which(non_founder==i)
                matrix.person[[3 * i - 2]] <- mutation_paternal
                matrix.person[[3 * i - 1]] <- mutation_maternal
                names(dimnames(matrix.person[[3 * i - 2]])) <- list(paste(arc.mutation[2*I-1,2],"_L1",sep=""),paste(arc.mutation[2*I-1,1],"_L1",sep=""))
                names(dimnames(matrix.person[[3 * i - 1]])) <- list(paste(arc.mutation[2*I,2],"_L1",sep=""),paste(arc.mutation[2*I,1],"_L1",sep=""))
                matrix.parent_to_child[[length(matrix.parent_to_child)+1]] <- parent_to_child
                names(dimnames(matrix.parent_to_child[[length(matrix.parent_to_child)]])) <- list(paste(arc.parent_to_child[6*I-3,2],"_L1",sep=""),
                paste(arc.parent_to_child[6*I-4,1],"_L1",sep=""),paste(arc.parent_to_child[6*I-5,1],"_L1",sep=""),paste(arc.parent_to_child[6*I-3,1],"_L1",sep=""))
                matrix.parent_to_child[[length(matrix.parent_to_child)+1]] <- parent_to_child
                names(dimnames(matrix.parent_to_child[[length(matrix.parent_to_child)]])) <- list(paste(arc.parent_to_child[6*I,2],"_L1",sep=""),
                paste(arc.parent_to_child[6*I-1,1],"_L1",sep=""),paste(arc.parent_to_child[6*I-2,1],"_L1",sep=""),paste(arc.parent_to_child[6*I,1],"_L1",sep=""))
              }
              for (i in 1:n){
                matrix.person[[3*i]] <- allele_to_genotype
                names(dimnames(matrix.person[[3*i]])) <- list(paste(node.person2[3 * i],"_L1",sep=""),paste(node.person2[3 * i - 2],"_L1",sep=""),paste(node.person2[3 * i - 1],"_L1",sep=""))
              }
              for(i in 1:nrow(arc.L1_to_L2)){
                matrix.linkage[[i]] <- matrix(c(0.5,0.5), ncol = 2, dimnames = list(NULL, Boolean))
              }
    
              matrix.list_L1 <- c(matrix.person,matrix.target,matrix.parent_to_child,matrix.linkage)
    
    
            }else if(k==linkage.loci[l,2]){
            
              names(dimnames(matrix.target[[1]])) <- list(paste(node.target[1],"_L2",sep=""),paste(node.person[3 * target_number - 2],"_L2",sep=""),"target")
              names(dimnames(matrix.target[[2]])) <- list(paste(node.target[2],"_L2",sep=""),paste(node.person[3 * target_number - 1],"_L2",sep=""),"target")
            
              for (i in founder){
                matrix.person[[3 * i - 2]] <- matrix(cpList2, ncol = length(aList2), dimnames = list(paste(node.person2[3*i-2],"_L2",sep=""), aList2))
                matrix.person[[3 * i - 1]] <- matrix(cpList2, ncol = length(aList2), dimnames = list(paste(node.person2[3*i-1],"_L2",sep=""), aList2))
              } 
              for (i in non_founder){
                I <- which(non_founder==i)
                matrix.person[[3 * i - 2]] <- mutation_paternal
                matrix.person[[3 * i - 1]] <- mutation_maternal      
                names(dimnames(matrix.person[[3 * i - 2]])) <- list(paste(arc.mutation[2*I-1,2],"_L2",sep=""),paste(arc.mutation[2*I-1,1],"_L2",sep=""))
                names(dimnames(matrix.person[[3 * i - 1]])) <- list(paste(arc.mutation[2*I,2],"_L2",sep=""),paste(arc.mutation[2*I,1],"_L2",sep=""))
                matrix.parent_to_child[[length(matrix.parent_to_child)+1]] <- parent_to_child
                names(dimnames(matrix.parent_to_child[[length(matrix.parent_to_child)]])) <- list(paste(arc.parent_to_child[6*I-3,2],"_L2",sep=""),
                paste(arc.parent_to_child[6*I-4,1],"_L2",sep=""),paste(arc.parent_to_child[6*I-5,1],"_L2",sep=""),paste(arc.parent_to_child[6*I-3,1],"_L2",sep=""))
                matrix.parent_to_child[[length(matrix.parent_to_child)+1]] <- parent_to_child
                names(dimnames(matrix.parent_to_child[[length(matrix.parent_to_child)]])) <- list(paste(arc.parent_to_child[6*I,2],"_L2",sep=""),
                paste(arc.parent_to_child[6*I-1,1],"_L2",sep=""),paste(arc.parent_to_child[6*I-2,1],"_L2",sep=""),paste(arc.parent_to_child[6*I,1],"_L2",sep=""))
              } 
              for (i in 1:n){
                matrix.person[[3*i]] <- allele_to_genotype
                names(dimnames(matrix.person[[3*i]])) <- list(paste(node.person2[3 * i],"_L2",sep=""),paste(node.person2[3 * i - 2],"_L2",sep=""),paste(node.person2[3 * i - 1],"_L2",sep=""))
              } 
              for(i in 1:nrow(arc.L1_to_L2)){
                matrix.linkage[[i]] <- L1_L2
                names(dimnames(matrix.linkage[[i]])) <- list(arc.L1_to_L2[i,2],arc.L1_to_L2[i,1])
              }
          
              matrix.list_L2 <- c(matrix.person,matrix.target,matrix.parent_to_child,matrix.linkage)
            }
          }
        
          matrix_list <- c(matrix.list_L1,matrix.list_L2,list(target_query))
          names(matrix_list) <- c(nodenames_L1,nodenames_L2,"target")
          dfit <- custom.fit(PES,dist = matrix_list)
  
          fitted.grain <-  as.grain(dfit)
          known_gtnodes <- c(nodenames_L1[3 * known_numbers],nodenames_L2[3 * known_numbers])
          known_gt_L1 <- apply(alleletype[,linkage.loci[l,1],DNAtype,m,drop=FALSE],3,function(x){return(paste(min(x),"-",max(x),sep=""))})
          known_gt_L2 <- apply(alleletype[,linkage.loci[l,2],DNAtype,m,drop=FALSE],3,function(x){return(paste(min(x),"-",max(x),sep=""))})
    
          net2 <- setEvidence(fitted.grain, nodes=known_gtnodes, states=c(known_gt_L1,known_gt_L2)) 
          result.part <- querygrain(net2,nodes=c("target"))
          result.LR.L[m,l] <- result.part$target["T"] / result.part$target["F"]
        }
      }
    }
    result.LR[[1]] <- result.LR.noL
    result.LR[[2]] <- result.LR.L
    return(result.LR)
  }

  Calculate <- function(Hypo.result){
    t <- proc.time()

    aList <- pList <- list()
    for (i in 1:(length(Hypo.result$locus_set))){
      aList[[i]] <-  Hypo.result$frequency_matrix[,1][!is.na(Hypo.result$frequency_matrix[,i+1])]
      pList[[i]] <-  Hypo.result$frequency_matrix[,i+1][!is.na(Hypo.result$frequency_matrix[,i+1])]
    }

    mutation_pat <- Hypo.result$mutation_matrix[,2:6]
    mutation_mat <- Hypo.result$mutation_matrix[,7:11]

    No.known.ref <- sum(Hypo.result$df[[1]]$affected)
    Known.ref <- Hypo.result$df[[1]]$DNAtype[which(Hypo.result$df[[1]]$DNAtype!=0)]

    alleletype <- array(NA, dim=c(2, length(Hypo.result$locus_set), No.known.ref, 1))
    for (i in 1:No.known.ref){
      alleletype[,,i,1] <- apply(t(Hypo.result$DNAtype_matrix[, (2 * i) : (2 * i + 1)]), c(1,2), function(x){return(as.numeric(x))})
    }

    df.original <- vector("list", length = 2)
    for (j in 1:2){
      df.original[[j]][[1]] <- pedigree.simple(Hypo.result$df[[j]])
    }

    LR.info <- list(aList, pList, mutation_pat, mutation_mat, alleletype, Hypo.result$linkage.loci, Hypo.result$RR)
    names(LR.info) <- c("aList", "pList", "mutation_pat", "mutation_mat", "alleletype", "linkage.loci", "RR")

    target_number <- c(NA, NA)
    check.result <- rep(NA,No.known.ref)

    for (i in Known.ref){
      df <- df.original
      for (j in 1:2){
        df[[j]][[1]]$DNAtype[which(df[[j]][[1]]$DNAtype==i)] <- 0
        df[[j]][[2]] <- pedigree.simple(df[[j]][[1]])
      }
      check.result[which(Known.ref==i)] <- check.pedigree(df[[1]][[2]],df[[2]][[2]])      
    }

    if(length(which(check.result==1))==1){      
      LR.step <- list()
      for (j in 1:2){
        target_number[j] <- which(check.result==1)
        if (is.na(Hypo.result$RR[1])){
          calc.info <- c(LR.info, list(df = df.original[[j]][[1]], target_number = target_number[j]))
          LR.step[[j]] <- PES_function(calc.info)
        }
        else{
          calc.info <- c(LR.info, list(df = df.original[[j]][[1]], target_number = target_number[j]))
          temp.result <- PES_function_linkage(calc.info)
          LR.step[[j]] <- list()
          LR.step[[j]][[1]] <- temp.result[[1]]
          LR.step[[j]][[2]] <- temp.result[[2]]
        }
      }

     }else{
       LR.step <- list()
       for (j in 1:2){ 
         if (is.na(Hypo.result$RR[1])){
           LR.step[[j]] <- array(NA, dim=c(No.known.ref, length(Hypo.result$locus_set)))
         }else{
           LR.step[[j]] <- list()
           LR.step[[j]][[1]] <- array(NA, dim=c(No.known.ref, length(Hypo.result$locus_set) - 2 * length(Hypo.result$RR)))
           LR.step[[j]][[2]] <- array(NA, dim=c(No.known.ref, length(Hypo.result$RR)))
         }
       }
       for (i in 1:No.known.ref){
         for (j in 1:2){
           target_number[j] <- which(df.original[[j]][[i]]$DNAtype == Known.ref[i])
           if (is.na(Hypo.result$RR[1])){
             calc.info <- c(LR.info, list(df = df.original[[j]][[i]], target_number = target_number[j]))
             LR.step[[j]][i,] <- PES_function(calc.info)
           }
           else{
             calc.info <- c(LR.info, list(df = df.original[[j]][[i]], target_number = target_number[j]))
             temp.result <- PES_function_linkage(calc.info)
             LR.step[[j]][[1]][i,] <- temp.result[[1]]
             LR.step[[j]][[2]][i,] <- temp.result[[2]]
           }
          if (i != No.known.ref){
            df.original[[j]][[i]]$DNAtype[target_number[j]] <- 0         
            df.original[[j]][[i + 1]] <- pedigree.simple(df.original[[j]][[i]])
          }
        }
      }
    }
   
    if (is.na(Hypo.result$RR[1])){
      LR <- apply(LR.step[[1]], 2, prod) / apply(LR.step[[2]], 2, prod)
    }
    else{
      LR <- list()
      LR[[1]] <- apply(LR.step[[1]][[1]], 2, prod) / apply(LR.step[[2]][[1]], 2, prod)
      LR[[2]] <- apply(LR.step[[1]][[2]], 2, prod) / apply(LR.step[[2]][[2]], 2, prod)
    }

    calc.time <- proc.time() - t
    return_information <- list(LR, calc.time)
    names(return_information) <- c("LR","calc.time")
    return(c(Hypo.result, return_information))
  }

  Pedigree.make_linkage_mutation <- function(pedmake.info){
    pedigree <- apply(pedmake.info$df[,c("pat","mat")], c(1,2), pedigree.numbering, pedmake.info$df[,"id"])
    al <- pedmake.info$aList
    pl <- pedmake.info$pList

    result <- array(0,dim=c(2,length(al),nrow(pedigree),pedmake.info$N))
    for (j in 1:pedmake.info$N){
      i <- 1
      while(is.element(0, result[1,1,,j])){
        if (result[1,1,i,j] == 0){
          if (pedigree[i,1] == 0){
            result[,,i,j] <- mapply(sample,al,2,prob=pl,replace=TRUE)
          }
          else if (pedigree[i,1] != 0 && result[1,1,pedigree[i,1],j] != 0 && result[1,1,pedigree[i,2],j] != 0){
            if (!is.na(pedmake.info$RR[1])){
              for(k in 1:length(al)){
                if(any(k==pedmake.info$linkage.loci[,1])){
                  K <- which(pedmake.info$linkage.loci[,1]==k)     
                  haplotype <- sample(c(1,2),2,replace=T)
                  result[1,k,i,j] <- result[haplotype[1],k,pedigree[i,1],j]
                  result[2,k,i,j] <- result[haplotype[2],k,pedigree[i,2],j]
                  result[1,pedmake.info$linkage.loci[K,2],i,j] <- result[sample(c(haplotype[1],3-haplotype[1]),1,prob=c(1-pedmake.info$RR[K],pedmake.info$RR[K])),pedmake.info$linkage.loci[K,2],pedigree[i,1],j]
                  result[2,pedmake.info$linkage.loci[K,2],i,j] <- result[sample(c(haplotype[2],3-haplotype[2]),1,prob=c(1-pedmake.info$RR[K],pedmake.info$RR[K])),pedmake.info$linkage.loci[K,2],pedigree[i,2],j]
                }else if(all(k!=pedmake.info$linkage.loci[,2])){
                  result[,k,i,j] <- apply(result[,k,pedigree[i,],j],2,sample,1)
                }         
              }
            }
            else{
              result[,,i,j] <- apply(result[,,pedigree[i,],j],c(3,2),sample,1)
            }
            for(k in 1:length(al)){
              result[1,k,i,j] <- result[1,k,i,j] + sample(c(-2,-1,0,1,2),1,prob=pedmake.info$mutation_pat[k,])
              result[2,k,i,j] <- result[2,k,i,j] + sample(c(-2,-1,0,1,2),1,prob=pedmake.info$mutation_mat[k,])
              result[which(min(al[[k]]) > result[,k,i,j]),k,i,j] <-  min(al[[k]])
              result[which(max(al[[k]]) < result[,k,i,j]),k,i,j] <-  max(al[[k]]) 
            }
          }
        }
        if (i == nrow(pedigree)){
          i <- 1
        }
        else{
          i <- i + 1
        }
      }
    }
    return(result)
  }

  Simulate <- function(Hypo.result){
    t <- proc.time()

    aList <- pList <- list()
    for (i in 1:(length(Hypo.result$locus_set))){
      aList[[i]] <-  Hypo.result$frequency_matrix[,1][!is.na(Hypo.result$frequency_matrix[,i+1])]
      pList[[i]] <-  Hypo.result$frequency_matrix[,i+1][!is.na(Hypo.result$frequency_matrix[,i+1])]
    }

    mutation_pat <- Hypo.result$mutation_matrix[,2:6]
    mutation_mat <- Hypo.result$mutation_matrix[,7:11]

    No.known.ref <- sum(Hypo.result$df[[1]]$affected)
    Known.ref <- Hypo.result$df[[1]]$DNAtype[which(Hypo.result$df[[1]]$DNAtype!=0)]

    alleletype <- list()
    for (k in 1:2){
      pedmake.info <- list(aList, pList, mutation_pat, mutation_mat, Hypo.result$linkage.loci, Hypo.result$RR, Hypo.result$N, Hypo.result$df[[k]])
      names(pedmake.info) <- c("aList", "pList", "mutation_pat", "mutation_mat", "linkage.loci", "RR", "N", "df")
      alleletype[[k]] <- Pedigree.make_linkage_mutation(pedmake.info)[,,1:sum(Hypo.result$df[[k]]$affected),]
    }

    df.original <- vector("list", length = 2)
    for (j in 1:2){
      df.original[[j]][[1]] <- pedigree.simple(Hypo.result$df[[j]])
    }

    LR.info <- list()
    for (k in 1:2){
      LR.info[[k]] <- list(aList, pList, mutation_pat, mutation_mat, alleletype[[k]], Hypo.result$linkage.loci, Hypo.result$RR)
      names(LR.info[[k]]) <- c("aList", "pList", "mutation_pat", "mutation_mat", "alleletype", "linkage.loci", "RR")
    }

    target_number <- c(NA, NA)
    check.result <- rep(NA,No.known.ref)

    for (i in Known.ref){
      df <- df.original
      for (j in 1:2){
        df[[j]][[1]]$DNAtype[which(df[[j]][[1]]$DNAtype==i)] <- 0
        df[[j]][[2]] <- pedigree.simple(df[[j]][[1]])
      }
      check.result[which(Known.ref==i)] <- check.pedigree(df[[1]][[2]],df[[2]][[2]])      
    }

    LR.step <- vector("list", length = 2) 
   
    if(length(which(check.result==1))==1){
      for (k in 1:2){
        for(j in 1:2){
          target_number[j] <- which(df.original[[j]][[1]]$DNAtype == which(check.result==1))
          if (is.na(Hypo.result$RR[1])){
            LR.step[[k]][[j]] <- array(NA, dim=c(Hypo.result$N, length(Hypo.result$locus_set), 1))
            calc.info <- c(LR.info[[k]], list(df = df.original[[j]][[1]], target_number = target_number[j]))
            LR.step[[k]][[j]][,,1] <- PES_function(calc.info)
          }else{
            LR.step[[k]][[j]] <- list()
            LR.step[[k]][[j]][[1]] <- array(NA, dim=c(Hypo.result$N, length(Hypo.result$locus_set) - 2 * length(Hypo.result$RR),1))
            LR.step[[k]][[j]][[2]] <- array(NA, dim=c(Hypo.result$N, length(Hypo.result$RR), 1))
            calc.info <- c(LR.info[[k]], list(df = df.original[[j]][[1]], target_number = target_number[j]))
            temp.result <- PES_function_linkage(calc.info)
            LR.step[[k]][[j]][[1]][,,1] <- temp.result[[1]]
            LR.step[[k]][[j]][[2]][,,1] <- temp.result[[2]]
          }
        }
      }      
    }else{          
      for (k in 1:2){        
        for(j in 1:2){
          if (is.na(Hypo.result$RR[1])){
            LR.step[[k]][[j]] <- array(NA, dim=c(Hypo.result$N, length(Hypo.result$locus_set), sum(Hypo.result$df[[j]]$affected)))  
          }else{
            LR.step[[k]][[j]] <- list()
            LR.step[[k]][[j]][[1]] <- array(NA, dim=c(Hypo.result$N, length(Hypo.result$locus_set) - 2 * length(Hypo.result$RR), sum(Hypo.result$df[[1]]$affected)))
            LR.step[[k]][[j]][[2]] <- array(NA, dim=c(Hypo.result$N, length(Hypo.result$RR), sum(Hypo.result$df[[j]]$affected)))
          }               
          for (i in 1:sum(Hypo.result$df[[j]]$affected)){
            target_number[j] <- which(df.original[[j]][[i]]$DNAtype == i)
            calc.info <- c(LR.info[[k]], list(df = df.original[[j]][[i]], target_number = target_number[j]))            
            if (is.na(Hypo.result$RR[1])){          
              LR.step[[k]][[j]][,,i] <- PES_function(calc.info)
            }else{
              temp.result <- PES_function_linkage(calc.info)
              LR.step[[k]][[j]][[1]][,,i] <- temp.result[[1]]
              LR.step[[k]][[j]][[2]][,,i] <- temp.result[[2]]
            }
            if (i != sum(Hypo.result$df[[j]]$affected)){
              df.del <- df.original[[j]][[i]]
              df.del$DNAtype[target_number[j]] <- 0
              df.original[[j]][[i + 1]] <- pedigree.simple(df.del)
            }
          }
        }
      }
    }

    if (is.na(Hypo.result$RR[1])){
      LR <- list()
      for (k in 1:2){
        LR[[k]] <- apply(LR.step[[k]][[1]], c(1,2), prod) / apply(LR.step[[k]][[2]], c(1,2), prod)
      }
    }else{
      LR <- vector("list", length = 2)
      for (k in 1:2){
        LR[[k]][[1]] <- apply(LR.step[[k]][[1]][[1]], c(1,2), prod) / apply(LR.step[[k]][[2]][[1]], c(1,2), prod)
        LR[[k]][[2]] <- apply(LR.step[[k]][[1]][[2]], c(1,2), prod) / apply(LR.step[[k]][[2]][[2]], c(1,2), prod)
      }
    }

    calc.time <- proc.time() - t
    return_information <- list(LR, calc.time)
    names(return_information) <- c("LR","calc.time")
    return(c(Hypo.result, return_information))
  }

  Report.common <- function(result){
    report_filepath.frequency <- tclvalue(result$filepath_frequency)
    report_filepath.mutation <- tclvalue(result$filepath_mutation)
    if (result$mode == "LR"){
      report_filepath.DNAtype <- tclvalue(result$filepath_DNAtype)
    }
    report_locus_set_name <- result$locus_set_name
    report_locus_set <- paste(result$locus_set, sep="", collapse=",")

    result_pedigree.array <- list()
    for (j in 1:2){
      result_pedigree.array[[j]] <- array("-", dim=c(nrow(result$df[[j]]), 3))
      for (i in 1:nrow(result$df[[j]])){
        result_pedigree.array[[j]][i, 1] <- as.character(result$df[[j]]$names)[i]
        if (result$df[[j]]$pat[i] != 0){
          result_pedigree.array[[j]][i, 2] <- as.character(result$df[[j]]$names)[result$df[[j]]$pat[i]]
        }
        if (result$df[[j]]$mat[i] != 0){
          result_pedigree.array[[j]][i, 3] <- as.character(result$df[[j]]$names)[result$df[[j]]$mat[i]]
        }
      }
    }
    report_pedigree <- paste("H1\nName,Father,Mother\n",
                         paste(apply(result_pedigree.array[[1]], 1, function(x){paste(x[1], x[2], x[3], sep=",")}), sep="", collapse="\n"),"\n",
                         "\nH2\nName,Father,Mother\n",
                         paste(apply(result_pedigree.array[[2]], 1, function(x){paste(x[1], x[2], x[3], sep=",")}), sep="", collapse="\n"), sep=""
                            )

    if (!is.na(result[["RR"]][1])){
      report.linkage.loci <- apply(result$linkage.loci,c(1,2),function(x){return(result[["locus_set"]][x])})
      report.linkage.allelenames <- paste("Locus 1,Locus 2,Recombination rate\n", paste(paste(apply(report.linkage.loci, 1, paste, collapse=","), result[["RR"]], sep=","), collapse="\n"), sep="", collapse="")
    }
    else{
      report.linkage.allelenames <- "No linkage"
    }
    if (result$mode == "LR"){
      report_common <- paste(c("======== Files ========\n",
                          "Allele frequencies :,",report_filepath.frequency,"\n",
                          "Mutation rates :,",report_filepath.mutation,"\n",
                          "Profiles :,",report_filepath.DNAtype,"\n",
                          "\n======== Locus set ========\n",
                          "Kit :,",report_locus_set_name,"\n",
                          "Locus :,",report_locus_set,"\n",
                          "\n======== Linkage ========\n",
                          report.linkage.allelenames,"\n",
                          "\n======== Hypothesis ========\n",
                          report_pedigree
                        ),sep="",collapse=""
                      )
    }
    else{
      report_common <- paste(c("======== Files ========\n",
                          "Allele frequencies :,",report_filepath.frequency,"\n",
                          "Mutation rates :,",report_filepath.mutation,"\n",
                          "\n======== Locus set ========\n",
                          "Kit :,",report_locus_set_name,"\n",
                          "Locus :,",report_locus_set,"\n",
                          "\n======== Linkage ========\n",
                          report.linkage.allelenames,"\n",
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
        }
        else{
          file.create(SaveFileName)
          write(report,file=SaveFileName)
        }
      }
      else{
        SaveFileName <- paste(tclvalue(pre_SaveFileName),".csv",sep="")
        if (file.exists(SaveFileName)){
          warningMessage <- paste(strsplit(SaveFileName,"/")[[1]][length(strsplit(SaveFileName,"/")[[1]])],"? already exists.\nDo you want to overwrite it?")
          overwrite <- tkmessageBox(message=warningMessage, type="okcancel", icon="warning")
          if (tclvalue(overwrite) == "ok"){
            write(report,file=SaveFileName)
          }
        }
        else{
          file.create(SaveFileName)
          write(report,file=SaveFileName)
        }
      }
    }
  }

  LR_Report <- function(LR_result){
    if (!is.na(LR_result$RR[1])){
      linkage.loci.vector <- as.vector(LR_result$linkage.loci)
      report_LRprod <- prod(LR_result$LR[[1]]) * prod(LR_result$LR[[2]])
      report_result <- paste(LR_result$locus_set[-linkage.loci.vector], ":,", LR_result$LR[[1]], sep="", collapse="\n")

      report_linkage.loci <- apply(LR_result$linkage.loci,c(1,2),function(x){return(LR_result$locus_set[x])})
      report_linkage <- paste(paste(report_linkage.loci[,1], " - ", report_linkage.loci[,2], ":,", LR_result$LR[[2]]), collapse="\n")
    }
    else{
      report_result <- paste(paste(LR_result$locus_set, ":,", LR_result$LR, sep=""), sep="", collapse="\n")
      report_LRprod <- prod(LR_result$LR)
      report_linkage <- ""
    }
    report_today <- format(LR_result$today, "%Y/%b/%d %X")
    report_calc.time <- paste(signif(LR_result$calc.time[3], digits=3), "sec", collapse="")

    report_LR.specific <- paste(c("\n\n======== Likelihood ratio ========\n",
                              "Overall LR :, ",report_LRprod,"\n\n",
                              report_result,"\n\n", 
                              report_linkage,"\n",
                              "\n======== Information ========\n",
                              "Version: ,", version, "\n",
                              "Date: ,",report_today,"\n",
                              "Computation time: ,",report_calc.time
                            ),sep="",collapse=""
                          )
    report <- paste(Report.common(LR_result), report_LR.specific, sep="", collapse="")

    pre_SaveFileName <- tkgetSaveFile(parent = window, filetypes = "{{CSV files} {.csv}}")
    Report.save(report, pre_SaveFileName)
  }

  simu_Report <- function(simu_result){
    simu_items <- c("","min","5th percentile","25th percentile","50th percentile","75th percentile","95th percentile","max")
    report_today <- format(simu_result$today, "%Y/%b/%d %X")
    report_calc.time <- paste(signif(simu_result$calc.time[3], digits=3), "sec", collapse="")

    LRprod <- list()
    for (k in 1:2){
      if (!is.na(simu_result$RR[1])){
        LRprod[[k]] <- sort(apply(simu_result$LR[[k]][[1]],1,prod) * apply(simu_result$LR[[k]][[2]],1,prod))
      }
      else{
        LRprod[[k]] <- sort(apply(simu_result$LR[[k]],1,prod))
      }
    }

    report_simu<- vector("list", length = 2)
    for (k in 1:2){
      report_simu[[k]][[1]] <- paste(simu_items[2], quantile(LRprod[[k]])["0%"], sep=":,")
      report_simu[[k]][[2]] <- paste(simu_items[3], LRprod[[k]][floor(length(LRprod[[k]]) / 20) + 1], sep=":,")
      report_simu[[k]][[3]] <- paste(simu_items[4], quantile(LRprod[[k]])["25%"], sep=":,")
      report_simu[[k]][[4]] <- paste(simu_items[5], quantile(LRprod[[k]])["50%"], sep=":,")
      report_simu[[k]][[5]] <- paste(simu_items[6], quantile(LRprod[[k]])["75%"], sep=":,")
      report_simu[[k]][[6]] <- paste(simu_items[7], LRprod[[k]][ceiling(length(LRprod[[k]]) * 19 / 20)], sep=":,")
      report_simu[[k]][[7]] <- paste(simu_items[8], quantile(LRprod[[k]])["100%"], sep=":,")
    }

    report_simu.specific <- paste(c("\n\n======== Likelihood ratio ========\n",
                              "H1 true","\n",
                              paste(report_simu[[1]], collapse="\n"),
                              "\n\n","H2 true","\n",
                              paste(report_simu[[2]], collapse="\n"),
                              "\n\n======== Parameters ========\n",
                              "Number of simulations :,",simu_result$N,"\n",
                              "\n======== Information ========\n",
                              "Version: ,", version, "\n",
                              "Date: ,",report_today,"\n",
                              "Computation time: ,",report_calc.time
                            ),sep="",collapse=""
                          )
    report <- paste(Report.common(simu_result), report_simu.specific, sep="", collapse="")

    pre_SaveFileName <- tkgetSaveFile(parent = window, filetypes = "{{CSV files} {.csv}}")
    Report.save(report, pre_SaveFileName)
  }

  alldata_Report <- function(simu_result){
    alldata <- list()
    if (!is.na(simu_result$RR[1])){
      linkage.loci.vector <- as.vector(simu_result$linkage.loci)
      alldata.linkage.loci <- apply(simu_result$linkage.loci, c(1,2), function(x){return(simu_result$locus_set[x])})
      alldata.allelenames <- paste(c(simu_result$locus_set[-linkage.loci.vector], paste(alldata.linkage.loci[,1], " - ", alldata.linkage.loci[,2])), sep="", collapse=",")
      alldata.LRprod <- list()
      for (k in 1:2){
        alldata.LRprod[[k]] <- apply(simu_result$LR[[k]][[1]], 1, prod) * apply(simu_result$LR[[k]][[2]], 1, prod)
        alldata[[k]] <- paste("", apply(simu_result$LR[[k]][[1]], 1 ,paste, sep="", collapse=","), apply(simu_result$LR[[k]][[2]], 1, paste, sep="", collapse=","), alldata.LRprod[[k]], sep=",", collapse="\n")
      }
    }
    else{
      alldata.allelenames <- paste(simu_result$locus_set, sep="", collapse=",")
      alldata.LRprod <- list()
      for (k in 1:2){
        alldata.LRprod[[k]] <- apply(simu_result$LR[[k]], 1, prod)
        alldata[[k]] <- paste("", apply(simu_result$LR[[k]], 1, paste, sep="", collapse=","), alldata.LRprod[[k]], sep=",", collapse="\n")
      }
    }

    alldata <- paste(c(",", alldata.allelenames, ",", "Overall LR",
                     "\nH1 true",
                     alldata[[1]],"\n",
                     "\nH2 true",
                     alldata[[2]],"\n"
                      ), sep="", collapse=""
                    )

    pre_SaveFileName <- tkgetSaveFile(parent = window, filetypes = "{{CSV files} {.csv}}")
    Report.save(alldata, pre_SaveFileName)
  }


  if (sum(as.character(tkfont.names()) == "our_font")){
    tkfont.configure("our_font",family = "Meiryo",size = "10")
  }
  else{
    tkfont.create("our_font",family="Meiryo",size="10")
  }

  if (sum(as.character(tkfont.names()) == "bold_font")){
    tkfont.configure("bold_font",family = "Meiryo",size = "10",weight="bold")
  }
  else{
    tkfont.create("bold_font",family="Meiryo",size="10",weight="bold")
  }


  window <- tktoplevel()
  tkwm.title(window, "KinBN")

  Tabs <- tk2notebook(window,tabs = c("Mode", "Files", "Linkage", "Hypothesis", "Result (case)", "Result (simulation)"))
  tkpack(Tabs, fill = "both", expand = 1)
  tab.Mode <- tk2notetab(Tabs, "Mode")
  tab.Files <- tk2notetab(Tabs, "Files")
  tab.Linkage <- tk2notetab(Tabs, "Linkage")
  tab.Hypo <- tk2notetab(Tabs, "Hypothesis")
  tab.LR <- tk2notetab(Tabs, "Result (case)")
  tab.Simu <- tk2notetab(Tabs, "Result (simulation)")

  File.result <- list()
  Hypo.result <- list()
  LR_result <- list()
  simu_result <- list()
  linkage_number <- 0
  
  already_pack.File <- F
  already_pack.Linkage_Hypo <- F
  already_pack.apply <- F
  already_pack.LR <- F
  already_pack.simu <- F
  frame_Files <- tkframe(tab.Files)
  frame_Linkage <- tkframe(tab.Linkage)
  frame_Hypo <- tkframe(tab.Hypo)
  frame_LR <- tkframe(tab.LR)
  tkpack(frame_LR)
  frame_LR_result <- tkframe(frame_LR)
  frame_LR_simu <- tkframe(frame_LR)
  frame_Simu <- tkframe(tab.Simu)

  frame_Mode <- tkframe(tab.Mode)
  frame_title <- tkframe(frame_Mode)
  frame_mode.select <- tkframe(frame_Mode)

  mode <- tclVar("1")
  button_LR <- tkbutton(frame_mode.select,font="our_font",text="Case analysis", width=18, 
                         command=function(){
                           if (!already_pack.File || tclvalue(tkmessageBox(message="Setting of Linkage and Hypothesis will be deleted. Do you want to continue?", type="okcancel", parent = window, icon="warning")) == "ok"){
                             if (already_pack.File){
                               tkdestroy(frame_Files)
                               tkdestroy(frame_Linkage)
                               tkdestroy(frame_Hypo)
                               tkdestroy(frame_LR_result)
                               tkdestroy(frame_LR_simu)
                               tkdestroy(frame_Simu)
                               already_pack.Linkage_Hypo <<- already_pack.apply <<- already_pack.LR <<- already_pack.simu <<- F
                             }
                             else{
                               already_pack.File <<- T
                             }
                             frame_Files <<- tkframe(tab.Files)
                             tkpack(frame_Files)
                             makeFiles("LR")
                             tk2notetab.select(Tabs, "Files")
                           }
                                           }
                         )
  button_simu <- tkbutton(frame_mode.select,font="our_font",text="Simulation", width=18, 
                         command=function(){
                           if (!already_pack.File || tclvalue(tkmessageBox(message="Setting of Linkage and Hypothesis will be deleted. Do you want to continue?", type="okcancel", parent = window, icon="warning")) == "ok"){
                             if (already_pack.File){
                               tkdestroy(frame_Files)
                               tkdestroy(frame_Linkage)
                               tkdestroy(frame_Hypo)
                               tkdestroy(frame_LR_result)
                               tkdestroy(frame_LR_simu)
                               tkdestroy(frame_Simu)
                               already_pack.Linkage_Hypo <<- already_pack.apply <<- already_pack.LR <<- already_pack.simu <<- F
                             }
                             else{
                               already_pack.File <<- T
                             }
                             frame_Files <<- tkframe(tab.Files)
                             tkpack(frame_Files)
                             makeFiles("simu")
                             tk2notetab.select(Tabs, "Files")
                           }
                                           }
                         )
  tkpack(frame_Mode)
  tkpack(frame_mode.select,anchor="w",padx="30",pady="20",expand="1")

  tkpack(button_LR,anchor="w", side="left", padx=20)
  tkpack(button_simu,anchor="w")
}

KinBN()