checkFile <- function(envProj, modeStateVar){
  freqFp <- get("freqFp", pos = envProj)
  mutationFp <- get("mutationFp", pos = envProj)
  gtFp <- get("gtFp", pos = envProj)
  RR <- get("RR", pos = envProj)
  linkage.loci <- get("linkage.loci", pos = envProj)
  linkage_number <- get("linkage_number", pos = envProj)
  locus_set_name <- get("locus_set_name", pos = envProj)

  linkState <- get("linkState", pos = envProj)
  dropState <- get("dropState", pos = envProj)
  freqState <- get("freqState", pos = envProj)
  linkStateVar <- tclVar(linkState)
  dropStateVar <- tclVar(dropState)
  freqStateVar <- tclVar(freqState)
  DropMode <- get("DropMode", pos = envProj)
  LinkMode <- get("LinkMode", pos = envProj)
  LinkModeVar <- tclVar(LinkMode)
  default_locus_set <- get("default_locus_set", pos = envProj)
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
  assign("default_locus_set", default_locus_set, envir = envProj)

  if (length(freqFp) == 0){
    tkmessageBox(message = "Please load allele frequency file.", icon="error")
  }else if(length(mutationFp) == 0){
    tkmessageBox(message = "Please load mutation rates file.", icon="error")
  }else if(modeStateVar == "1" && length(gtFp) == 0){
    tkmessageBox(message = "Please load profiles file.", icon="error")
  }else{
    frequency_matrix <- read.csv(freqFp)
    locus_set <- names(frequency_matrix)[-1]
    mutation_matrix <- read.csv(mutationFp, stringsAsFactors=F)
    locus_set.mutation <- mutation_matrix[,1]
    if (modeStateVar == "1"){
      DNAtype_matrix <- read.csv(gtFp, stringsAsFactors=F)
      locus_set.DNAtype <- DNAtype_matrix[,1]
      allele.DNAtype <- DNAtype_matrix[,-1]
    }else{
      locus_set.DNAtype <- ""
      DNAtype_matrix <- ""
    }

    if(modeStateVar == "1" && any(is.na(allele.DNAtype))){
      tclvalue(dropStateVar) <-  "normal"
      assign("dropState", tclvalue(dropStateVar), envir = envProj)
      assign("DropMode", 2, envir = envProj)
    }

    if (ncol(mutation_matrix) != 11){
      tkmessageBox(message = "Error: incorrect file format of Mutation rate", icon="error")
    }else if (!identical(locus_set, locus_set.mutation)){
      tkmessageBox(message = "Error: incorrect locus set1", icon="error")
    }else if (modeStateVar == "1" && prod(sapply(locus_set, function(x){return(sum(locus_set.DNAtype == x) == 1)})) == 0 ){
      tkmessageBox(message = "Error: incorrect locus set2", icon="error")
    }else{
      is.default <- sapply(default_locus_set, function(x){if (setequal(x[[1]], locus_set) && length(x[[1]]) == length(locus_set) ) return(1) else return(0)})
      if (length(which(is.default == 1)) > 0 && tclvalue(LinkModeVar) == "1"){
        locus_set_name <- names(default_locus_set)[which(is.default == 1)]
        default.sort <- sapply(default_locus_set[[locus_set_name]][[1]], function(x){return(which(x == locus_set))})
        linkage.loci <- apply(default_locus_set[[locus_set_name]][[2]], c(1,2), function(x){return(default.sort[x])})
        RR <- default_locus_set[[locus_set_name]][[3]]
        linkage_number <- nrow(linkage.loci)
      }else if(length(which(is.default == 1)) > 0 && tclvalue(LinkModeVar) == "2"){
        locus_set_name <- names(default_locus_set)[which(is.default == 1)]
        default.sort <- sapply(default_locus_set[[locus_set_name]][[1]], function(x){return(which(x == locus_set))})
      }else if(length(which(is.default == 1)) > 0 && tclvalue(LinkModeVar) == "3"){
        locus_set_name <- names(default_locus_set)[which(is.default == 1)]
        default.sort <- sapply(default_locus_set[[locus_set_name]][[1]], function(x){return(which(x == locus_set))})
      }else{
        locus_set_name <- "-"
      }

      if (modeStateVar == "1"){
        DNAtype.sort <- sapply(locus_set, function(x){return(which(x == locus_set.DNAtype))})
        DNAtype_matrix <- DNAtype_matrix[DNAtype.sort, ]
      }

      if (modeStateVar == "2" || length(locus_set) == length(locus_set.DNAtype) || tclvalue(tkmessageBox(message="Loci without allele frequencies and mutation rates will be excluded from calculation.\nDo you want to continue?", type="okcancel", icon="warning")) == "ok"){

        assign("frequency_matrix", frequency_matrix, envir = envProj)
        assign("mutation_matrix", mutation_matrix, envir = envProj)
        assign("DNAtype_matrix", DNAtype_matrix, envir = envProj)
        assign("locus_set_name", locus_set_name, envir = envProj)
        assign("locus_set", locus_set, envir = envProj)
        assign("linkage.loci", linkage.loci, envir = envProj)
        assign("RR", RR, envir = envProj)
        assign("linkage_number", linkage_number, envir = envProj)
        assign("LinkMode", tclvalue(LinkModeVar), envir = envProj)
        tclvalue(linkStateVar) <-  "normal"
        tclvalue(freqStateVar) <-  "normal"
        assign("linkState", tclvalue(linkStateVar), envir = envProj)
        assign("freqState", tclvalue(freqStateVar), envir = envProj)
      }
    }
  }
}
