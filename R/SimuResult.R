SimuResult <- function(envProj, envGUI){

  tf <- tktoplevel()
  tkwm.title(tf, "Simulation result")
  outframe_simu <- tkframe(tf)
  frame_simu.result <- tkframe(outframe_simu)
  frame_simu.button <- tkframe(outframe_simu)

  return_information.simu <- list()
  return_information.simu <- get("return_information.simu", pos = envProj)

  LR.simu <- return_information.simu[[1]]

  softVer <- get("softVer", pos = envProj)
  locus_set <- get("locus_set", pos = envProj)
  linkage.loci <- get("linkage.loci", pos = envProj)
  RR <- get("RR", pos = envProj)
  OverallLR <- get("LRprod", pos = envProj)

  LRprod <- list()
  for (k in 1:2){
    if (!is.na(RR[1])){
      LRprod[[k]] <- sort(apply(LR.simu[[k]][[1]],1,prod) * apply(LR.simu[[k]][[2]],1,prod))
    }else{
      LRprod[[k]] <- sort(apply(LR.simu[[k]],1,prod))
    }
  }

  if(length(which(LRprod[[2]]!=0)) < 2 ){
    xmax <- max(c(density(log10(LRprod[[1]]))[["x"]]), log10(OverallLR))
    xmin <- min(c(density(log10(LRprod[[1]]))[["x"]]), log10(OverallLR))
    xwidth <- xmax - xmin
    xmax <- xmax + xwidth / 15
    xmin <- xmin - xwidth / 15
    ymax <- max(c(density(log10(LRprod[[1]]))[["y"]]) * 1.1)
    graph.draw <- function(){
      par(mar=c(4,6,3,9))
      plot(density(log10(LRprod[[1]])), col=rgb(1,0,0), lwd=2, las=1, cex.lab=1.2, xlim=c(xmin,xmax), ylim=c(0,ymax), xlab=expression(log[10](LR)), main="")
      par(new=T)
      abline(v=log10(OverallLR), col=rgb(0,0,0), lwd=3, lty=2)
      par(xpd=T)
      legend(par()$usr[2], par()$usr[4],col=c(rgb(1,0,0),rgb(0,0.5,1),rgb(0,0,0)),lty=c(1,1,2),lwd=c(3,3,3),c("H1 true","H2 true", "Case"))
    }
  }else if(length(which(!is.infinite(LRprod[[1]]))) < 2 ){
    xmax <- max(c(density(log10(LRprod[[2]]))[["x"]]), log10(OverallLR))
    xmin <- min(c(density(log10(LRprod[[2]]))[["x"]]), log10(OverallLR))
    xwidth <- xmax - xmin
    xmax <- xmax + xwidth / 15
    xmin <- xmin - xwidth / 15
    ymax <- max(c(density(log10(LRprod[[2]]))[["y"]]) * 1.1)
    graph.draw <- function(){
      par(mar=c(4,6,3,9))
      plot(density(log10(LRprod[[2]])), col=rgb(0,0.5,1), lwd=2, las=1, cex.lab=1.2, xlim=c(xmin,xmax), ylim=c(0,ymax), xlab=expression(log[10](LR)), main="")
      par(new=T)
      abline(v=log10(OverallLR), col=rgb(0,0,0), lwd=3, lty=2)
      par(xpd=T)
      legend(par()$usr[2], par()$usr[4],col=c(rgb(1,0,0),rgb(0,0.5,1),rgb(0,0,0)),lty=c(1,1,2),lwd=c(3,3,3),c("H1 true","H2 true", "Case"))
    }
  }else{
    xmax <- max(c(density(log10(LRprod[[1]]))[["x"]], density(log10(LRprod[[2]]))[["x"]]))
    xmin <- min(c(density(log10(LRprod[[1]]))[["x"]], density(log10(LRprod[[2]]))[["x"]]))
    xwidth <- xmax - xmin
    xmax <- xmax + xwidth / 15
    xmin <- xmin - xwidth / 15
    ymax <- max(c(density(log10(LRprod[[1]]))[["y"]], density(log10(LRprod[[2]]))[["y"]])) * 1.1
    graph.draw <- function(){
      par(mar=c(4,6,3,9))
      plot(density(log10(LRprod[[1]])), col=rgb(1,0,0), lwd=2, xlim=c(xmin,xmax), ylim=c(0,ymax), xlab="", ylab="", ann=F, axes=F, main="")
      par(new=T)
      plot(density(log10(LRprod[[2]])), col=rgb(0,0.5,1), lwd=2, las=1, cex.lab=1.2, xlim=c(xmin,xmax), ylim=c(0,ymax), xlab=expression(log[10](LR)), main="")
      abline(v=log10(OverallLR), col=rgb(0,0,0), lwd=3, lty=2)
      par(xpd=T)
      legend(par()$usr[2], par()$usr[4],col=c(rgb(1,0,0),rgb(0,0.5,1),rgb(0,0,0)),lty=c(1,1,2),lwd=c(3,3,3),c("H1 true","H2 true", "Case"))
    }
  }

  simu_graph <- tkrplot(frame_simu.result, function() graph.draw(), hscale=2, vscale=1.2)

  button_simu.report <- tkbutton(frame_simu.button,text=" Report ", cursor="hand2",
                                 command=function() simu_Report(envProj))
  button_simu.download <- tkbutton(frame_simu.button,text=" All data ", cursor="hand2",
                                   command=function() alldata_Report(envProj))
  button_simu.graph <- tkbutton(frame_simu.button,text=" Open R graphics device ", cursor="hand2",
                                 command=function() graph.draw())

  tkpack(simu_graph)
  tkpack(button_simu.report,anchor="w",side="left",padx="10",pady="10")
  tkpack(button_simu.download,anchor="w",side="left",padx="10",pady="10")
  tkpack(button_simu.graph,anchor="w",side="left",padx="10",pady="10")
  tkpack(frame_simu.result,anchor="w")
  tkpack(frame_simu.button,anchor="w")
  tkpack(outframe_simu)
  assign("outframe_simu", outframe_simu, envir = envGUI)
}
