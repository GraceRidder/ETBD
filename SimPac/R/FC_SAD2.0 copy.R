#SAD generators
library(sads)

MakeSAD = function(matlist, sadmarg, Jmax){

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

#  f <- function(S) {
#    Jmax * (S / (100 + S))
#  }


  J <- c()
 for ( i in 1:length(Jmax)){
  J <- append(J, Jmax[i] * (S[i] / (100 + S[i])))
 }

 # J <- f(SV)


#Brent solver
alphaBrentsolver=function(S,
                          J,
                          bounds=c(0.01,100) #bounds of search interval for alpha
)
{
  #defining the alpha function
  fn=function(alpha) alpha*log(1+(J/alpha))-S

  #looking for more appropriate bounds that have opposite signs of functional value
  if (sign(fn(bounds[1]))+sign(fn(bounds[2]))!=0){
    repeat{
      bounds[1]=bounds[1]/10
      bounds[2]=bounds[2]*10
      if(sign(fn(bounds[1]))+sign(fn(bounds[2]))==0) break
    }
  }

  #solver
  root=uniroot(fn, bounds)

  results=list(alpha=root$root, #alpha
               bounds=bounds #final bounds
  )
  return(results)
}



alpha <- c()
for (i in 1:length(matlist)){
  alpha <- append(alpha, alphaBrentsolver(S = S[i], J = J[i], bounds=c(0.01,100))$alpha)
}








#iterative stochastic simulator of logseries given alpha and

happyls_sads=function(S,
                      J,
                      alpha,
                      reltol #relative tolerance of resulting J
)
{
  repeat{
    newabund=rls(n = S,N = J,alpha = alpha)
    if(sum(newabund) < J+(J*reltol) && sum(newabund) > J-(J*reltol)) break
  }
  return(newabund)
}


newabund <- matlist
for (i in 1:length(matlist)) {
  newabund[[i]] = happyls_sads(S[i], J[i], alpha[i], reltol)
  newabund[[i]] <- sort(newabund[[i]], decreasing = T)
}



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






