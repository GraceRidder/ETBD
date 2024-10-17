#SAD generators
library(sads)

MakeNOSAD = function(matlist, Jmax){
  
  
  matnames <- unmatrixlist(matlist)
  matab <- unlist(matlist)
  S <- length(matnames)
  nmatab <- -matab
  
  SV <- c()
  for (i in 1:length(matlist)) {
    SV <- append(SV, length(matlist[[i]]))
  }
  
  S <- SV
  
  J <- c()
  for ( i in 1:length(Jmax)){
    J <- append(J, Jmax[i] * (S[i] / (100 + S[i])))
  }
  
  Jax <- J
  Sax <- S
  newb <- list()
  
MED <- J/S


for (i in 1:length(Jmax)) {
  J <- Jax[i]
  S <- Sax[i]
  MED1 <- MED[i]
  newabund <- rep(MED1, S)
  newb[[i]] <- newabund
  
}


  mat <- matlist
  for (i in 1:length(matlist)){
    for (k in 1:length(matlist[[i]])){
      mat[[i]][k] <- newb[[i]][k]
    }
  }
  

  return(mat)
  
}

