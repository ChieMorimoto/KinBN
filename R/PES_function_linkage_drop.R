PES_function_linkage_drop <- function(calc.info.all, DropMode, FixedLR, PrD, FreqMode, min.freq, locus.number){
  DNAtype <- calc.info.all$df[,"DNAtype"]
  alleletype <- calc.info.all$alleletype
  target_number <- calc.info.all$target_number
  pa_mutation <- calc.info.all$mutation_pat
  ma_mutation <- calc.info.all$mutation_mat
  al <- calc.info.all$aList
  pl <- calc.info.all$pList
  range <- 0

  linkage.loci <- calc.info.all$linkage.loci
  RR <- calc.info.all$RR
  linkage.loci.vector <- as.vector(t(linkage.loci))
  locusdrop.loci <- integer(0)

  drop.loci <-which(apply(alleletype, 2,function(x){any(is_blank(x))})==TRUE)
  locusdrop.loci <- which(apply(apply(alleletype, c(3,2),function(x){all(is_blank(x))}),2,any)==TRUE)
  drop.loci <- drop.loci[which(is.na(match(drop.loci,locusdrop.loci))==TRUE)]

  result.LR.locus <- matrix(1, 1, length(al))

  calc.info <- calc.info.all

  if(length(locusdrop.loci) > 0){
    linkage.loci.locusdrop <- (which(linkage.loci.vector %in% locusdrop.loci==TRUE) %/% 2) + (which(linkage.loci.vector %in% locusdrop.loci==TRUE) %% 2)
  }

  if(DropMode == "0" && any(drop.loci %in% linkage.loci.vector)){
    linkage.loci.D <- unique((which(linkage.loci.vector %in% drop.loci==TRUE) %/% 2) + (which(linkage.loci.vector %in% drop.loci==TRUE) %% 2))
    locus.number <- locus.number[!(locus.number %in% linkage.loci[linkage.loci.D,])]
    result.LR.locus <- PES_function_linkage(calc.info, FreqMode, min.freq, locus.number)
  }else if(DropMode == "0" && !any(drop.loci %in% linkage.loci.vector)){
    result.LR.locus <- PES_function_linkage(calc.info, FreqMode, min.freq, locus.number)
  }else if(length(linkage.loci.vector)==0){
    result.LR.locus <- PES_function_drop(calc.info, DropMode, FixedLR, PrD, FreqMode, min.freq, locus.number)
  }else if(!any(drop.loci %in% linkage.loci.vector)){
    locus.number.noL <- locus.number[!(locus.number %in% linkage.loci.vector)]
    result.LR.locus <- PES_function_drop(calc.info, DropMode, FixedLR, PrD, FreqMode, min.freq, locus.number.noL)
    linkloci <- as.vector(calc.info$linkage.loci[,1])
    A <- PES_function_linkage(calc.info, FreqMode, min.freq, linkloci)
    result.LR.locus[1,calc.info$linkage.loci[,1]] <- A[linkloci]
  }else if(any(drop.loci %in% linkage.loci.vector)){
    locus.number.noL <- locus.number[!(locus.number %in% linkage.loci.vector)]
    result.LR.locus <- PES_function_drop(calc.info, DropMode, FixedLR, PrD, FreqMode, min.freq, locus.number.noL)

    linkage.loci.D <- unique((which(linkage.loci.vector %in% drop.loci==TRUE) %/% 2) + (which(linkage.loci.vector %in% drop.loci==TRUE) %% 2))
    if(length(linkage.loci.D) != nrow(linkage.loci)){
      locus.number.L <- locus.number[!(locus.number %in% linkage.loci[linkage.loci.D,])]
      result.LR.locus[1,calc.info$linkage.loci[-linkage.loci.D,1]] <- PES_function_linkage(calc.info, FreqMode, min.freq, locus.number.L)[calc.info$linkage.loci[-linkage.loci.D,1]]
    }

    pedigree <- apply(calc.info$df[,c("pat","mat")], c(1,2), pedigree.numbering, calc.info$df[,"id"])
    range <- 0
    n <- dim(pedigree)[1]
    founder <- which(pedigree[,1] == 0)
    non_founder <- which(pedigree[,1] != 0)
    known_numbers <- which(calc.info$df$DNAtype != 0)

    linkdrop <-  c(1:(nrow(linkage.loci)))
    if(length(locusdrop.loci) > 0){
      linkdrop <- which(!(c(1:(nrow(linkage.loci))) %in% linkage.loci.locusdrop)==TRUE)
    }

    if (length(non_founder) >= 1){
      for (l in linkdrop){
        nodekind <- c("pg","mg","gt")
        node.person <- unlist(lapply(paste("n",1:n,sep=""),function(x){return(paste(x,nodekind,sep=""))}))
        node.target <- c(paste("p",target_number,c("pg","mg"),sep=""))
        node.mutation <- unlist(lapply(paste("o",non_founder,sep=""),function(x){return(paste(x,c("pg","mg"),sep=""))}))
        node.linkage <- unlist(lapply(paste("o",non_founder,sep=""),function(x){return(paste(x,c("pg_fp","mg_fp"),sep=""))}))

        nodenames <- c(node.person,node.target,node.mutation,node.linkage)
        nodenames_L1 <- c(paste(nodenames,"_L1",sep=""))
        nodenames_L2 <- c(paste(nodenames,"_L2",sep=""))

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

        matrix.list_L1 <- list()
        matrix.list_L2 <- list()

        arc.drop.L1 <- integer(0)
        arc.drop.L2 <- integer(0)

        known_gtnodes <- integer(0)

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

          Drop_person <- which(apply(alleletype[,k,,], 2,function(x){any(is_blank(x))})==TRUE)

          if(k==linkage.loci[l,1] && length(Drop_person)==0){

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
            known_gtnodes_L1 <- c(nodenames_L1[3 * known_numbers])
            known_gtnodes <- c(known_gtnodes,nodenames_L1[3 * known_numbers])


          }else if(k==linkage.loci[l,1] && length(Drop_person)>0){

            drop_person_L1 <- which(calc.info$df$DNAtype == Drop_person)
            node.drop <- unlist(lapply(paste("d",drop_person_L1,sep=""),function(x){return(paste(x,"gt",sep=""))}))

            nodenames1 <- c(node.person,node.target,node.mutation,node.linkage,node.drop)
            nodenames2 <- c(node.person,node.target,node.mutation,node.linkage)
            nodenames_L1 <- c(paste(nodenames1,"_L1",sep=""))
            nodenames_L2pre <- c(paste(nodenames2,"_L2",sep=""))


            allele1 <- (1:pol^2 - 1) %/% pol
            allele2 <- (1:pol^2 - 1) %% pol

            homo <- which(allele1 == allele2)
            gList.drop <- gList
            gList.drop[homo] <- unlist(lapply(aList2,function(x){return(paste(x,"-",NA,sep=""))}))
            gList2.drop <- c(gList.drop[-del],"NA-NA")

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

            matrix.drop <- list()

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
            for (i in 1:length(drop_person_L1)){
              matrix.drop[[i]] <- genotype_to_dropgenotype
              arc.drop.L1 <- matrix(nrow=length(drop_person_L1),ncol=2)
              arc.drop.L1[,1] <- as.vector(paste(sapply(drop_person_L1,function(x){return(node.person[3 * x])}),"_L1",sep=""))
              arc.drop.L1[,2] <- paste(node.drop,"_L1",sep="")
              names(dimnames(matrix.drop[[i]])) <- list(arc.drop.L1[i,2],arc.drop.L1[i,1])
            }
            matrix.list_L1 <- c(matrix.person,matrix.target,matrix.parent_to_child,matrix.linkage,matrix.drop)
            known_gtnodes_L1 <- c(nodenames_L1[3 * known_numbers])
            known_gtnodes_L1[which(known_numbers==drop_person_L1)] <- c(tail(nodenames_L1,n=length(drop_person_L1)))
            known_gtnodes <- c(known_gtnodes,known_gtnodes_L1)

          }else if(k==linkage.loci[l,2] && length(Drop_person)==0){

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
            known_gtnodes_L2 <- c(nodenames_L2[3 * known_numbers])
            known_gtnodes <- c(known_gtnodes,nodenames_L2[3 * known_numbers])

          }else if(k==linkage.loci[l,2] && length(Drop_person)>0){
            drop_person_L2 <- which(calc.info$df$DNAtype == Drop_person)
            node.drop <- unlist(lapply(paste("d",drop_person_L2,sep=""),function(x){return(paste(x,"gt",sep=""))}))

            nodenames1 <- c(node.person,node.target,node.mutation,node.linkage)
            nodenames2 <- c(node.person,node.target,node.mutation,node.linkage,node.drop)
            nodenames_L1pre <- c(paste(nodenames1,"_L1",sep=""))
            nodenames_L2 <- c(paste(nodenames2,"_L2",sep=""))

            allele1 <- (1:pol^2 - 1) %/% pol
            allele2 <- (1:pol^2 - 1) %% pol

            homo <- which(allele1 == allele2)
            gList.drop <- gList
            gList.drop[homo] <- unlist(lapply(aList2,function(x){return(paste(x,"-",NA,sep=""))}))
            gList2.drop <- c(gList.drop[-del],"NA-NA")

            if(DropMode == "3"){
              genotype_to_dropgenotype <- matrix(0,nrow=length(gList2.drop),ncol=length(gList2))
              homo_gt <- which(gList2 %in% gList[which(allele1 == allele2)])
              hetero_gt <- c(1:length(gList2))[-homo_gt]
              for (i in homo_gt){
                genotype_to_dropgenotype[i,i] <- 1
              }
              for (i in hetero_gt){
                genotype_to_dropgenotype[i,i] <- 1 - as.numeric(PrD)
                genotype_to_dropgenotype[homo_gt[(allele1[-del][i]+1)],i] <- as.numeric(PrD)*0.5
                genotype_to_dropgenotype[homo_gt[(allele2[-del][i]+1)],i] <- as.numeric(PrD)*0.5
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

            matrix.drop <- list()

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
            for (i in 1:length(drop_person_L2)){
              matrix.drop[[i]] <- genotype_to_dropgenotype
              arc.drop.L2 <- matrix(nrow=length(drop_person_L2),ncol=2)
              arc.drop.L2[,1] <- as.vector(paste(sapply(drop_person_L2,function(x){return(node.person[3 * x])}),"_L2",sep=""))
              arc.drop.L2[,2] <- paste(node.drop,"_L2",sep="")

              names(dimnames(matrix.drop[[i]])) <- list(arc.drop.L2[i,2],arc.drop.L2[i,1])
            }
            matrix.list_L2 <- c(matrix.person,matrix.target,matrix.parent_to_child,matrix.linkage,matrix.drop)
            known_gtnodes_L2 <- c(nodenames_L2[3 * known_numbers])
            known_gtnodes_L2[which(known_numbers==drop_person_L2)] <- c(tail(nodenames_L2,n=length(drop_person_L2)))
            known_gtnodes <- c(known_gtnodes,known_gtnodes_L2)
          }
        } #k
        PES <- empty.graph(c(nodenames_L1,nodenames_L2,"target"))
        arc.set <- rbind(arc.putative_or_true, arc.set_L1, arc.set_L2, arc.L1_to_L2,arc.drop.L1,arc.drop.L2)
        arcs(PES) <- arc.set

        matrix_list <- c(matrix.list_L1,matrix.list_L2,list(target_query))
        names(matrix_list) <- c(nodenames_L1,nodenames_L2,"target")
        dfit <- custom.fit(PES,dist = matrix_list)


        known_gt_L1 <- apply(alleletype[,linkage.loci[l,1],DNAtype,1,drop=FALSE],3,allele_to_dropgenotype)
        known_gt_L2 <- apply(alleletype[,linkage.loci[l,2],DNAtype,1,drop=FALSE],3,allele_to_dropgenotype)

        fitted.grain <-  as.grain(dfit)

        net2 <- setEvidence(fitted.grain, nodes=known_gtnodes, states=c(known_gt_L1,known_gt_L2))
        result.part <- querygrain(net2,nodes=c("target"))
        result.LR.locus[1,calc.info$linkage.loci[l,1]]  <- result.part$target["T"] / result.part$target["F"]
      }
    }
  }
  return(result.LR.locus)
}
