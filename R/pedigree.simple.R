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
