calcLR <- function(envProj, envGUI){

  labelframe_LR <- get("labelframe_LR", pos = envGUI)
  simuframe_LR <- get("simuframe_LR", pos = envGUI)
  finResults <- get("finResults", pos = envProj)
  if(finResults){
    inputOk <- tclvalue(tkmessageBox(message="Previous results will be deleted. Do you want to continue?",
                                     type="okcancel", icon="warning"))
    if(inputOk == "ok"){
      tkdestroy(labelframe_LR)
      tkdestroy(simuframe_LR)
    }
  }

  t <- proc.time()
  today <- Sys.time()
  df <- list()
  df[[1]] <- get("df_Hypo1", pos = envProj)
  df[[2]] <- get("df_Hypo2", pos = envProj)
  DNAtype_matrix <- get("DNAtype_matrix", pos = envProj)
  frequency_matrix <- get("frequency_matrix", pos = envProj)
  mutation_matrix <- get("mutation_matrix", pos = envProj)
  locus_set_name <- get("locus_set_name", pos = envProj)
  locus_set <- get("locus_set", pos = envProj)
  linkage.loci <- get("linkage.loci", pos = envProj)
  RR <- get("RR", pos = envProj)
  DropMode <- get("DropMode", pos = envProj)
  FixedLR <- get("FixedLR", pos = envProj)
  PrD <- get("PrD", pos = envProj)
  allele.DNAtype <- DNAtype_matrix[,-1]

  FreqMode <- get("FreqMode", pos = envProj)
  min.freq <- get("min.freq", pos = envProj)

  aList <- pList <- list()
  for (i in 1:(length(locus_set))){
    aList[[i]] <-  frequency_matrix[,1][!is.na(frequency_matrix[,i+1])]
    pList[[i]] <-  frequency_matrix[,i+1][!is.na(frequency_matrix[,i+1])]
  }

  mutation_pat <- mutation_matrix[,2:6]
  mutation_mat <- mutation_matrix[,7:11]

  No.known.ref <- sum(df[[1]]$affected)
  Known.ref <- df[[1]]$DNAtype[which(df[[1]]$DNAtype!=0)]

  alleletype <- array(NA, dim=c(2, length(locus_set), No.known.ref, 1))
  for (i in 1:No.known.ref){
    alleletype[,,i,1] <- apply(t(DNAtype_matrix[, (2 * i) : (2 * i + 1)]), c(1,2), function(x){return(as.numeric(x))})
  }

  drop.loci <-which(apply(alleletype, 2,function(x){any(is_blank(x))})==TRUE)
  locusdrop.loci <- which(apply(apply(alleletype, c(3,2),function(x){all(is_blank(x))}),2,any)==TRUE)
  drop.loci <- drop.loci[which(is.na(match(drop.loci,locusdrop.loci))==TRUE)]

  df.original <- vector("list", length = 2)
  for (j in 1:2){
    df.original[[j]][[1]] <- pedigree.simple(df[[j]])
  }

  LR.info <- list(aList, pList, mutation_pat, mutation_mat, alleletype, linkage.loci, RR)
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
      if (is.na(RR[1]) && any(is.na(allele.DNAtype))){
        calc.info <- c(LR.info, list(df = df.original[[j]][[1]], target_number = target_number[j]))
        LR.step[[j]] <- PES_function_drop(calc.info, DropMode, FixedLR, PrD, FreqMode,min.freq, c(1:length(locus_set)))
      }else if (is.na(RR[1])){
        calc.info <- c(LR.info, list(df = df.original[[j]][[1]], target_number = target_number[j]))
        LR.step[[j]] <- PES_function(calc.info,FreqMode,min.freq, c(1:length(locus_set)))

      }else if(any(is.na(allele.DNAtype))){
        calc.info <- c(LR.info, list(df = df.original[[j]][[1]], target_number = target_number[j]))
        LR.step[[j]] <- PES_function_linkage_drop(calc.info, DropMode, FixedLR, PrD,FreqMode,min.freq, c(1:length(locus_set)))
      }else{
        calc.info <- c(LR.info, list(df = df.original[[j]][[1]], target_number = target_number[j]))
        LR.step[[j]] <- PES_function_linkage(calc.info, FreqMode, min.freq, c(1:length(locus_set)))
      }
    }
  }else{
    LR.step <- list()
    for (j in 1:2){
      LR.step[[j]] <- array(NA, dim=c(No.known.ref, length(locus_set)))
    }
    for (i in 1:No.known.ref){
      for (j in 1:2){
        target_number[j] <- which(df.original[[j]][[i]]$DNAtype == Known.ref[i])
        if(is.na(RR[1]) && any(is.na(allele.DNAtype))){
          calc.info <- c(LR.info, list(df = df.original[[j]][[i]], target_number = target_number[j]))
          LR.step[[j]][i,] <- PES_function_drop(calc.info, DropMode, FixedLR, PrD, FreqMode, min.freq, c(1:length(locus_set)))
        }else if(is.na(RR[1])){
          calc.info <- c(LR.info, list(df = df.original[[j]][[i]], target_number = target_number[j]))
          LR.step[[j]][i,] <- PES_function(calc.info, FreqMode, min.freq, c(1:length(locus_set)))
        }else if(any(is.na(allele.DNAtype))){
          calc.info <- c(LR.info, list(df = df.original[[j]][[i]], target_number = target_number[j]))
          LR.step[[j]][i,] <- PES_function_linkage_drop(calc.info, DropMode, FixedLR, PrD,FreqMode,min.freq, c(1:length(locus_set)))
        }else{
          calc.info <- c(LR.info, list(df = df.original[[j]][[i]], target_number = target_number[j]))
          LR.step[[j]][i,] <- PES_function_linkage(calc.info, FreqMode, min.freq, c(1:length(locus_set)))
        }
        if (i != No.known.ref){
          df.original[[j]][[i]]$DNAtype[target_number[j]] <- 0
          df.original[[j]][[i + 1]] <- pedigree.simple(df.original[[j]][[i]])
        }
      }
    }
  }
  if (is.na(RR[1])){
    LR <- apply(LR.step[[1]], 2, prod) / apply(LR.step[[2]], 2, prod)
    if(DropMode == "0" && length(drop.loci) > 0){
      LR[drop.loci] <- FixedLR
    }
  }else{
    LR <- apply(LR.step[[1]], 2, prod) / apply(LR.step[[2]], 2, prod)
    if(DropMode == "0" && length(drop.loci) > 0){
      LR[drop.loci] <- FixedLR
    }
  }
  calc.time <- proc.time() - t
  return_information <- list(LR, calc.time, today)
  names(return_information) <- c("LR","calc.time", "today")
  assign("return_information", return_information, envir = envProj)
  makeLR(envProj, envGUI)
}
