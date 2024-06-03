PES_function_linkage <- function(calc.info, FreqMode, min.freq, locus.number){
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
  linkage.loci.vector <- as.vector(t(linkage.loci))

  result.LR.locus <- matrix(1, 1, length(al))
  if(length(al) > length(as.vector(linkage.loci))){
    locus.number.noL <- locus.number[!(locus.number %in% linkage.loci.vector)]
    result.LR.locus <- PES_function(calc.info, FreqMode, min.freq, locus.number.noL)
  }

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

    for (l in 1:nrow(linkage.loci)){
      for (k in linkage.loci[l,]){
        aList_number <- sort(unique(as.vector(alleletype[,k,,1])))
        uni_aList <- unique(c(al[[k]],aList_number))
        aList <- sort(unique(as.vector(sapply(c(-range:range),function(x){return(x + uni_aList)}))))

        cp_convert <- function(x){
          if(FreqMode == "0"){
            if (any(al[[k]] == x)){
              return((pl[[k]][which(al[[k]] == x)] + 1) / (sum(pl[[k]]) + length(aList)))
            }else{
              return(1 / (sum(pl[[k]]) + length(aList)))
            }
          }else if(FreqMode == "1"){
            if (!any(al[[k]] == x)){
              return(5 / sum(pl[[k]]))
            }else if (pl[[k]][which(al[[k]] == x)] < 5){
              return(5 / sum(pl[[k]]))
            }else{
              return((pl[[k]][which(al[[k]] == x)] + 1) / (sum(pl[[k]]) + 1))
            }
          }else{
            if (!any(al[[k]] == x)){
              return(min.freq)
            }else{
              return(pl[[k]][which(al[[k]] == x)]) / (sum(pl[[k]]))
            }
          }
        }

        cpList <- sapply(aList,cp_convert)
        cpList <- (cpList/sum(cpList))


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
        mutation_paternal <- apply(pre_mutation_paternal,1,function(x){return(x / sum(x))})
        dimnames(mutation_paternal) <- list(aList2,aList2)

        pre_mutation_maternal <- matrix(0,nrow=pol,ncol=pol)
        for (i in 1:length(aList2)){
          for (g in -2:2){
            if (any(aList2 == aList2[i] + g)){
              pre_mutation_maternal[i,which(aList2 == aList2[i] + g)] <- ma_mutation[k,(g+3)]
            }
          }
        }
        mutation_maternal <- apply(pre_mutation_maternal,1,function(x){return(x / sum(x))})
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
      known_gt_L1 <- apply(alleletype[,linkage.loci[l,1],DNAtype,1,drop=FALSE],3,function(x){return(paste(min(x),"-",max(x),sep=""))})
      known_gt_L2 <- apply(alleletype[,linkage.loci[l,2],DNAtype,1,drop=FALSE],3,function(x){return(paste(min(x),"-",max(x),sep=""))})

      net2 <- setEvidence(fitted.grain, nodes=known_gtnodes, states=c(known_gt_L1,known_gt_L2))
      result.part <- querygrain(net2,nodes=c("target"))
      result.LR.locus[,linkage.loci[l,1]] <- result.part$target["T"] / result.part$target["F"]
    }
  }
  return(result.LR.locus)
}
