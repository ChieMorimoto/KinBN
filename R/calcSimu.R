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


calcSimu <- function(envProj, envGUI, LR.N, modeStateVar){

  labelframe_LR <- get("labelframe_LR", pos = envGUI)
  simuframe_LR <- get("simuframe_LR", pos = envGUI)
  outframe_simu <- get("outframe_simu", pos = envGUI)

  finResults <- get("finResults", pos = envProj)

  if(modeStateVar == "2" && finResults){
    inputOk <- tclvalue(tkmessageBox(message="Previous results will be deleted. Do you want to continue?",
                                     type="okcancel", icon="warning"))
    if(inputOk == "ok"){
      tkdestroy(labelframe_LR)
      tkdestroy(simuframe_LR)
      tkdestroy(outframe_simu)

    }
  }

  if (is.na(as.integer(LR.N))){
    tkmessageBox(message = "Please enter an integer of >2 to Number of simulation.", icon="error")
    stop()
  }
  else if (as.integer(LR.N) < 2){
    tkmessageBox(message = "Please enter an integer of >2 to Number of simulation.", icon="error")
    stop()
  }

  t.simu <- proc.time()
  today.simu <- Sys.time()
  df <- list()
  df[[1]] <- get("df_Hypo1", pos = envProj)
  df[[2]] <- get("df_Hypo2", pos = envProj)

  frequency_matrix <- get("frequency_matrix", pos = envProj)
  mutation_matrix <- get("mutation_matrix", pos = envProj)
  locus_set_name <- get("locus_set_name", pos = envProj)
  locus_set <- get("locus_set", pos = envProj)
  linkage.loci <- get("linkage.loci", pos = envProj)
  linkage.loci.vector <- as.vector(linkage.loci)
  RR <- get("RR", pos = envProj)

  FreqMode <- get("FreqMode", pos = envProj)
  min.freq <- get("min.freq", pos = envProj)


  t.simu <- proc.time()
  n.simu <- as.numeric(LR.N)
  pb <- tkProgressBar(title = "Simulation", min = 0, max = n.simu, width = 300)

  aList <- pList <- list()
  for (i in 1:(length(locus_set))){
    aList[[i]] <-  frequency_matrix[,1][!is.na(frequency_matrix[,i+1])]
    pList[[i]] <-  frequency_matrix[,i+1][!is.na(frequency_matrix[,i+1])]
  }

  mutation_pat <- mutation_matrix[,2:6]
  mutation_mat <- mutation_matrix[,7:11]

  No.known.ref <- sum(df[[1]]$affected)
  Known.ref <- df[[1]]$DNAtype[which(df[[1]]$DNAtype!=0)]

  alleletype <- list()
  for (k in 1:2){
    pedmake.info <- list(aList, pList, mutation_pat, mutation_mat, linkage.loci, RR, n.simu, df[[k]])
    names(pedmake.info) <- c("aList", "pList", "mutation_pat", "mutation_mat", "linkage.loci", "RR", "N", "df")
    alleletype[[k]] <- Pedigree.make_linkage_mutation(pedmake.info)[,,1:sum(df[[k]]$affected),]
  }

  df.original <- vector("list", length = 2)
  for (j in 1:2){
    df.original[[j]][[1]] <- pedigree.simple(df[[j]])
  }

  LR.info <- list()
  for (k in 1:2){
    LR.info[[k]] <- list(aList, pList, mutation_pat, mutation_mat, alleletype[[k]], linkage.loci, RR)
    names(LR.info[[k]]) <- c("aList", "pList", "mutation_pat", "mutation_mat", "alleletype", "linkage.loci", "RR")
  }

  target_number <- c(NA, NA)
  check.result <- rep(NA,No.known.ref)

  for (i in Known.ref){
    df.step <- df.original
    for (j in 1:2){
      df.step[[j]][[1]]$DNAtype[which(df.step[[j]][[1]]$DNAtype==i)] <- 0
      df.step[[j]][[2]] <- pedigree.simple(df.step[[j]][[1]])
    }
    check.result[which(Known.ref==i)] <- check.pedigree(df.step[[1]][[2]],df.step[[2]][[2]])
  }

  LR.step <- vector("list", length = 2)

  if(length(which(check.result==1))==1){
    for (k in 1:2){
      for(j in 1:2){
        LR.step[[k]][[j]] <- array(1, dim=c(n.simu, length(locus_set), 1))
      }
    }

  }else{
    for (k in 1:2){
      for(j in 1:2){
        LR.step[[k]][[j]] <- array(1, dim=c(n.simu, length(locus_set), sum(df[[j]]$affected)))
      }
    }
  }

  for(m in 1:n.simu){
    LR.info.each <- list()
    for (k in 1:2){
      LR.info.each[[k]] <- list(aList, pList, mutation_pat, mutation_mat, alleletype[[k]][,,,m,drop=F], linkage.loci, RR)
      names(LR.info.each[[k]]) <- c("aList", "pList", "mutation_pat", "mutation_mat", "alleletype", "linkage.loci", "RR")
    }

    if(length(which(check.result==1))==1){
      for (k in 1:2){
        for(j in 1:2){
          target_number[j] <- which(df.original[[j]][[1]]$DNAtype == which(check.result==1))
          if (is.na(RR[1])){
            calc.info <- c(LR.info.each[[k]], list(df = df.original[[j]][[1]], target_number = target_number[j]))
            LR.step[[k]][[j]][m,,1] <- PES_function(calc.info, FreqMode, min.freq, c(1:length(locus_set)))
          }else{
            calc.info <- c(LR.info.each[[k]], list(df = df.original[[j]][[1]], target_number = target_number[j]))
            LR.step[[k]][[j]][m,,1] <- PES_function_linkage(calc.info, FreqMode, min.freq, c(1:length(locus_set)))
          }
        }
      }

    }else{
      for (k in 1:2){
        for(j in 1:2){
          for (i in 1:sum(df[[j]]$affected)){
            target_number[j] <- which(df.original[[j]][[i]]$DNAtype == i)
            calc.info <- c(LR.info.each[[k]], list(df = df.original[[j]][[i]], target_number = target_number[j]))
            if (is.na(RR[1])){
              LR.step[[k]][[j]][m,,i] <- PES_function(calc.info, FreqMode, min.freq, c(1:length(locus_set)))
            }else{
              LR.step[[k]][[j]][m,,i] <- PES_function_linkage(calc.info, FreqMode, min.freq, c(1:length(locus_set)))
            }
            if (i != sum(df[[j]]$affected)){
              df.del <- df.original[[j]][[i]]
              df.del$DNAtype[target_number[j]] <- 0
              df.original[[j]][[i + 1]] <- pedigree.simple(df.del)
            }
          }
        }
      }
    }
    setTkProgressBar(pb,m, label=paste(round(m/n.simu*100, 0),"% done"))
  }
  close(pb)

  if (is.na(RR[1])){
    LR <- list()
    for (k in 1:2){
      LR[[k]] <- apply(LR.step[[k]][[1]], c(1,2), prod) / apply(LR.step[[k]][[2]], c(1,2), prod)
    }
  }else{
    LR <- vector("list", length = 2)

    for (k in 1:2){
      LR[[k]][[1]] <- apply(LR.step[[k]][[1]][,-linkage.loci.vector,], c(1,2), prod) / apply(LR.step[[k]][[2]][,-linkage.loci.vector,], c(1,2), prod)
      LR[[k]][[2]] <- matrix(1,n.simu,length(RR))
      for(i in 1:length(RR)){
        LR[[k]][[2]][,i] <- (apply(as.matrix(LR.step[[k]][[1]][,linkage.loci[i,1],]),1,prod) *
                               apply(as.matrix(LR.step[[k]][[1]][,linkage.loci[i,2],]),1,prod) ) /
          (apply(as.matrix(LR.step[[k]][[2]][,linkage.loci[i,1],]),1,prod) *
             apply(as.matrix(LR.step[[k]][[2]][,linkage.loci[i,2],]),1,prod))
      }
    }
  }

  calc.time.simu <- proc.time() - t.simu

  tkmessageBox(message = "Simulation calculations are complete.", icon = "info")

  return_information.simu <- list(LR, calc.time.simu, today.simu)
  names(return_information.simu) <- c("LR","calc.time.simu", "today.simu")
  assign("return_information.simu", return_information.simu, envir = envProj)
  assign("n.simu", n.simu, envir = envProj)
  makeSimu(envProj, envGUI, modeStateVar)
}
