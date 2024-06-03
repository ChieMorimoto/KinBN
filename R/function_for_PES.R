cp_convert <- function(x){
  if (any(al[[k]] == x)){
    return((pl[[k]][which(al[[k]] == x)] + 1) / (sum(pl[[k]]) + length(aList)))
  }
  else{
    return(1 / (sum(pl[[k]]) + length(aList)))
  }
}

is_blank <- function(x) {is.na(x) | x == ""}

allele_to_dropgenotype <- function(x){
  if(is.na(x[1])){
    x[1] <- x[2]
    x[2] <- NA
    return(paste(x[1],"-",x[2],sep=""))
  }else if(is.na(x[2])){
    return(paste(x[1],"-",x[2],sep=""))
  }
  return(paste(min(x),"-",max(x),sep=""))
}



