PES_function_drop <- function(calc.info, DropMode, FixedLR, PrD, FreqMode, min.freq, locus.number){
  pedigree <- apply(calc.info$df[,c("pat","mat")], c(1,2), pedigree.numbering, calc.info$df[,"id"])
  DNAtype <- calc.info$df[,"DNAtype"]
  alleletype <- calc.info$alleletype
  target_number <- calc.info$target_number
  pa_mutation <- calc.info$mutation_pat
  ma_mutation <- calc.info$mutation_mat
  al <- calc.info$aList
  pl <- calc.info$pList
  range <- 0

  result.LR <- matrix(1,1,length(al))

  drop.allloci <-which(apply(alleletype, 2,function(x){any(is_blank(x))})==TRUE)
  locusdrop.loci <- which(apply(apply(alleletype, c(3,2),function(x){all(is_blank(x))}),2,any)==TRUE)
  drop.loci <- drop.allloci[which(is.na(match(drop.allloci,locusdrop.loci))==TRUE)]

  n <- dim(pedigree)[1]
  founder <- which(pedigree[,1] == 0)
  non_founder <- which(pedigree[,1] != 0)
  known_numbers <- which(calc.info$df$DNAtype != 0)

  if(DropMode == "0" || length(drop.loci)==0){
    if(length(locusdrop.loci)>0){
      locus.number <- locus.number[!(locus.number %in% locusdrop.loci)]
      result.LR <- PES_function(calc.info, FreqMode, min.freq, locus.number)
    }else{
      result.LR <- PES_function(calc.info, FreqMode, min.freq, locus.number)
    }
  }else if(length(drop.loci)>0) {
    locus.number <- locus.number[!(locus.number %in% drop.loci)]
    result.LR <- PES_function(calc.info, FreqMode, min.freq, locus.number)
    for (k in drop.loci){
      Drop_person <- which(apply(alleletype[,k,,], 2,function(x){any(is_blank(x))})==TRUE)
      drop_person <- which(calc.info$df$DNAtype == Drop_person)

      if (length(non_founder) >= 1){
        nodekind <- c("pg","mg","gt")
        node.person <- unlist(lapply(paste("n",1:n,sep=""),function(x){return(paste(x,nodekind,sep=""))}))
        node.target <- c(paste("p",target_number,c("pg","mg"),sep=""),"target")
        node.mutation <- unlist(lapply(paste("o",non_founder,sep=""),function(x){return(paste(x,c("pg","mg"),sep=""))}))
        node.drop <- unlist(lapply(paste("d",drop_person,sep=""),function(x){return(paste(x,"gt",sep=""))}))

        nodenames <- c(node.person,node.target,node.mutation,node.drop)

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

        arc.drop <- matrix(nrow=length(drop_person),ncol=2)
        arc.drop[,1] <- as.vector(sapply(drop_person,function(x){return(node.person[3 * x])}))
        arc.drop[,2] <- node.drop

        arc.set <- rbind(arc.putative_or_true, arc.allele_to_genotype, arc.parent_to_child, arc.mutation,arc.drop)
        arcs(PES) <- arc.set


        aList_number <- sort(unique(as.vector(alleletype[,k,,1])))
        uni_aList <- unique(c(al[[k]],aList_number))
        range_aList <- sort(unique(as.vector(sapply(c(-range:range),function(x){return(x + uni_aList)}))))
        aList <- range_aList[which(min(uni_aList) <= range_aList & range_aList <= max(uni_aList))]
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
        allele1 <- (1:pol^2 - 1) %/% pol
        allele2 <- (1:pol^2 - 1) %% pol
        del <- which(allele1 > allele2)
        homo <- which(allele1 == allele2)
        gList.drop <- gList
        gList2 <- gList[-del]
        gList.drop[homo] <- unlist(lapply(aList2,function(x){return(paste(x,"-",NA,sep=""))}))
        gList2.drop <- c(gList.drop[-del],"NA-NA")

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

        genotype_to_dropgenotype <- matrix(0,nrow=length(gList2.drop),ncol=length(gList2))

        if(DropMode == "3"){
          genotype_to_dropgenotype <- matrix(0,nrow=length(gList2.drop),ncol=length(gList2))
          homo_gt <- which(gList2 %in% gList[which(allele1 == allele2)])
          hetero_gt <- c(1:length(gList2))[-homo_gt]
          Pr_D <- as.numeric(PrD)
          for (i in homo_gt){
            genotype_to_dropgenotype[i,i] <- 1-Pr_D^2
            genotype_to_dropgenotype[length(gList2.drop),i] <- Pr_D^2
          }
          for (i in hetero_gt){
            genotype_to_dropgenotype[i,i] <- (1 - Pr_D)^2
            genotype_to_dropgenotype[homo_gt[(allele1[-del][i]+1)],i] <- Pr_D*(1-Pr_D)
            genotype_to_dropgenotype[homo_gt[(allele2[-del][i]+1)],i] <- Pr_D*(1-Pr_D)
            genotype_to_dropgenotype[length(gList2.drop),i] <- Pr_D^2
          }
        }else if(DropMode == "2"){
          genotype_to_dropgenotype <- matrix(0,nrow=length(gList2.drop),ncol=length(gList2))
          homo_gt <- which(gList2 %in% gList[which(allele1 == allele2)])
          hetero_gt <- c(1:length(gList2))[-homo_gt]
          for (i in homo_gt){
            genotype_to_dropgenotype[i,i] <- 1/4
            genotype_to_dropgenotype[length(gList2.drop),i] <- 3/4
          }
          for (i in hetero_gt){
            genotype_to_dropgenotype[i,i] <- 1/4
            genotype_to_dropgenotype[homo_gt[(allele1[-del][i]+1)],i] <- 1/4
            genotype_to_dropgenotype[homo_gt[(allele2[-del][i]+1)],i] <- 1/4
            genotype_to_dropgenotype[length(gList2.drop),i] <- 1/4
          }
        }
        dimnames(genotype_to_dropgenotype) <- list(gList2.drop,gList2)

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
        matrix.mutation <- list()
        matrix.drop <- list()

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

        for (i in 1:length(drop_person)){
          matrix.drop[[i]] <- genotype_to_dropgenotype
          names(dimnames(matrix.drop[[i]])) <- list(arc.drop[i,2],arc.drop[i,1])
        }

        matrix_list <- c(matrix.person,matrix.target,matrix.mutation,matrix.drop)

        names(matrix_list) <- nodenames
        dfit <- custom.fit(PES,dist = matrix_list)
        fitted.grain <-  as.grain(dfit)

        known_gt <- apply(alleletype[,k,DNAtype,1,drop=FALSE],3,allele_to_dropgenotype)

        known_gtnodes <- c(nodenames[3 * known_numbers])
        known_gtnodes[which(known_numbers==drop_person)] <- c(tail(nodenames,n=length(drop_person)))
        net2 <- setEvidence(fitted.grain, nodes=known_gtnodes, states=known_gt)
        result.part <- querygrain(net2,nodes=c("target"))
        LR.locus <- result.part$target["T"] / result.part$target["F"]
        result.LR[1,k] <- LR.locus
      }
    }
  }
  return(result.LR)
}
