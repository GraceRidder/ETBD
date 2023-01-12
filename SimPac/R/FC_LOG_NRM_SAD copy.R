#SAD generators
library(sads)

MakeLogNormalSAD = function(matlist, sadmarg, Jmax){

  matnames <- unmatrixlist(matlist)
  matab <- unlist(matlist)
  reltol <- sadmarg
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

  for (i in 1:length(Jmax)){
    J <- Jax[i]
    S <- Sax[i]

####  log-normal
mm <- 0 # mean
## so I just start with a mean of zero and increase it until I get close..
#and then repeat the random generation until I'm within a 10% margin of J
newabund<- rlnorm(S, meanlog = mm, sdlog = 1)
sum(newabund)

if (sum(newabund) < J){
  repeat{
    mm <- mm + .01
    newabund<- rlnorm(S, meanlog = mm, sdlog = 1)
    if(sum(newabund) < J+(J*reltol) && sum(newabund) > J-(J*reltol)) break
  }
}
newb[[i]] <- newabund


}

  newabund<- newb

#attribute new sizes
### if species have the same abundance they are randomly placed in different ranks (within the original rank zone)
rankab <- matlist
speciesrank <- matlist
for (i in 1:length(matlist)){
  speciesrank[[i]] =rank(-matlist[[i]], ties.method = "random")
  rankab[[i]]<- newabund[[i]][speciesrank[[i]]]
}

mat <- matlist
for (i in 1:length(matlist)){
  for (k in 1:length(matlist[[i]])){
    mat[[i]][k] <- rankab[[i]][k]
  }
}


return(mat)

}





