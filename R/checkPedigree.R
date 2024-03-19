checkPedigree <- function(j){
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
    tkmessageBox(message = paste("Hypothesis", j, ": Error: incorrect setting of Sex", collapse="", sep=""), parent = window, icon="error")
  }else if (is.element(-1, df$pat) || is.element(-1, df$mat)){
    tkmessageBox(message = paste("Hypothesis", j, ": Error: incorrect setting of Father or Mother", collapse="", sep=""), parent = window, icon="error")
  }else{
    err <- try(pedigree(id = df$id, sex= df$sex, dadid = df$pat, momid = df$mat), silent = T)
    err.fc <- try(familycheck(famid=rep(1,nrow(df)), id = df$id, father.id = df$pat, mother.id = df$mat), silent = T)
    if (match(class(err),"try-error") || match(class(err.fc), "try-error")) {
      tkmessageBox(message = paste("Hypothesis", j, ": Error: incorrect setting of pedigree ", collapse="", sep=""), parent = window, icon="error")
    }else{
      return(df)
    }
  }
}

Hypo.pedigree <- function(df){
  if (is.data.frame(df)){
    window.ped <- tktoplevel()
    tkwm.title(window.ped, "View pedigree tree")
    unrelated <- which(df$pat == 0)[is.element(which(df$pat == 0), c(df$pat, df$mat)) == F]
    frame_ped <- tkframe(window.ped)
    tkpack(frame_ped)
    if (length(unrelated) != nrow(df)){
      labelframe_ped.related <- tk2labelframe(frame_ped, text="Pedigree tree", padding=c(20,20), labelanchor="nw", relief="groove", borderwidth="2")
      if (length(unrelated) != 0){
        labelframe_ped.related <- tk2labelframe(frame_ped, text="Pedigree tree", padding=c(20,20), labelanchor="nw", relief="groove", borderwidth="2")
        pedobj <- pedigree(id = df[-unrelated,]$id, sex= df[-unrelated,]$sex, dadid = df[-unrelated,]$pat, momid = df[-unrelated,]$mat, affected = df[-unrelated,]$affected)
        pedigree_graph <- tkRplot(labelframe_ped.related,function(){plot(pedobj, id=as.character(df[-unrelated,]$names))})
      }
      else{
        labelframe_ped.related <- tk2labelframe(frame_ped, text="Pedigree tree", padding=c(20,20), labelanchor="nw", relief="groove", borderwidth="2")
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
