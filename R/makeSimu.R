makeSimu <- function(envProj, envGUI, modeStateVar){
  return_information.simu <- list()
  return_information.simu <- get("return_information.simu", pos = envProj)

  LR.simu <- return_information.simu[[1]]

  softVer <- get("softVer", pos = envProj)
  locus_set <- get("locus_set", pos = envProj)
  linkage.loci <- get("linkage.loci", pos = envProj)
  RR <- get("RR", pos = envProj)

  LRprod <- list()
  for (k in 1:2){
    if (!is.na(RR[1])){
      LRprod[[k]] <- sort(apply(LR.simu[[k]][[1]],1,prod) * apply(LR.simu[[k]][[2]],1,prod))
    }else{
      LRprod[[k]] <- sort(apply(LR.simu[[k]],1,prod))
    }
  }
  if(modeStateVar == "1"){

    button_LR.simuResult <- get("button_LR.simuResult", pos = envGUI)
    tkdestroy(button_LR.simuResult)
    simuframe_LR <- get("simuframe_LR", pos = envGUI)
    button_LR.simuResult <- tkbutton(simuframe_LR, text="Result", cursor="hand2", width="10", state = "normal",
                                     command = function() SimuResult(envProj, envGUI))
    tkpack(button_LR.simuResult, side="left",anchor="n", padx=10, pady=10)
    tkpack(simuframe_LR, anchor="w", padx=10, pady=10)
    button_LR.simuResult <- assign("button_LR.simuResult", button_LR.simuResult, envir = envGUI)
    simuframe_LR <- assign("simuframe_LR", simuframe_LR, envir = envGUI)
  }else{
    tabs <- get("tabs", pos = envGUI)
    tabResults <- get("tabResults", pos = envGUI)
    labelframe_LR <- get("labelframe_LR", pos=envGUI)
    simuframe_LR <- get("simuframe_LR", pos=envGUI)
    tkdestroy(labelframe_LR)
    tkdestroy(simuframe_LR)

    outframe_simu <- tkframe(tabResults)
    frame_simu.result <- tkframe(outframe_simu)
    frame_simu.button <- tkframe(outframe_simu)

    if(length(which(LRprod[[2]]!=0)) < 2 ){
      xmax <- max(c(density(log10(LRprod[[1]]))[["x"]]))
      xmin <- min(c(density(log10(LRprod[[1]]))[["x"]]))
      xwidth <- xmax - xmin
      xmax <- xmax + xwidth / 15
      xmin <- xmin - xwidth / 15
      ymax <- max(c(density(log10(LRprod[[1]]))[["y"]]) * 1.1)
      graph.draw <- function(){
        par(mar=c(4,6,3,9))
        plot(density(log10(LRprod[[1]])), col=rgb(1,0,0), lwd=2, las=1, cex.lab=1.2, xlim=c(xmin,xmax), ylim=c(0,ymax), xlab=expression(log[10](LR)), main="")
        par(xpd=T)
        legend(par()$usr[2], par()$usr[4],col=c(rgb(1,0,0),rgb(0,0.5,1)),lty=c(1,1),lwd=c(3,3),c("H1 true","H2 true"))
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
        par(xpd=T)
        legend(par()$usr[2], par()$usr[4],col=c(rgb(1,0,0),rgb(0,0.5,1)),lty=c(1,1),lwd=c(3,3),c("H1 true","H2 true"))
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
    tk2notetab.select(tabs, "Results")
    assign("outframe_simu", outframe_simu, envir = envGUI)
    assign("finResults", TRUE, envir = envProj)
  }
}
