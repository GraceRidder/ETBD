###Geometric solver

#helper function for solving the geometric SAD
######################################################
geometricBrentsolver = function(S,
                                #number of species
                                J,
                                #total nuber of individuals
                                Nmin = 1,
                                #abundance of least abundant species
                                bounds = c(1.000001, 10),
                                #bounds of search interval for r
                                ...)
{
  #function to be solved
  fn = function(r)
    Nmin * ((1 - r ^ S) / (1 - r)) - J

  #looking for more appropriate bounds
  if (fn(bounds[1]) > 0) {
    repeat {
      bounds[1] = 1 + ((bounds[1] - 1) / 10)
      if (fn(bounds[1]) < 0)
        break
    }
  }

  #iterative solver using
  root = uniroot(fn, bounds, ...)

  #result object
  results = list(
    r = root$root,
    #r
    Nvect = Nmin * root$root ^ ((S:1) - 1),
    #vector of species abundances
    predict_J = (Nmin * (1 - root$root ^ S) / (1 - root$root)),
    #J prediction based on parameters (to check numerical errors)
    bounds = bounds #final bounds
  )

  #result
  return(results)

}





MakeGEO <- function(matlist, Jmax) {
  SV <- c()
  for (i in 1:length(matlist)) {
    SV <- append(SV, length(matlist[[i]]))
  }

  S <- SV


  J <- c()
  for (i in 1:length(matlist)) {
    J <- append(J, Jmax[i] * (S[i] / (100 + S[i])))
  }


  newabund <- matlist
  for (i in 1:length(matlist)) {
    if (S[i] > 5){
    newabund[[i]] = geometricBrentsolver(S = S[i],
                                         J = J[i],
                                         bounds = c(.01, 100))$Nvect
    }
  }

  #attribute new sizes
  ### if species have the same abundance they are randomly placed in different ranks (within the original rank zone)
  rankab <- matlist
  speciesrank <- matlist
  for (i in 1:length(matlist)) {
    speciesrank[[i]] = rank(-matlist[[i]], ties.method = "random")
    rankab[[i]] <- newabund[[i]][speciesrank[[i]]]
  }

  mat <- matlist
  for (i in 1:length(matlist)) {
    for (k in 1:length(matlist[[i]])) {
      mat[[i]][k] <- rankab[[i]][k]
    }
  }


  return(mat)

}




