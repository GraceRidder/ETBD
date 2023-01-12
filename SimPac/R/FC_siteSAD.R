## If sites are biomes, each site can have a SAD and Jmax

library(sads)


MakeSAD = function(matlist, sadmarg, Jmax){
  newmat <- matlist
  for (o in 1:length(matlist)){
  matnames <- row.names(matlist[[o]])
  reltol <-  sadmarg
  matab <- matlist[[o]][1:length(matnames)]
  S <- length(matnames)
  nmatab <- -matab
  f <- function(S) {
    Jmax[o]*(S/(100+S))
  }

  J <-f(S)

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

  alpha=alphaBrentsolver(S,J)$alpha


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


  newabund =happyls_sads(S, J, alpha, reltol)

  # S = 50
  # J = 300
  #
  # ####  log-normal
  # mm <- 0 # mean
  # ## so I just start with a mean of zero and increase it until I get close..
  # #and then repeat the random generation until I'm within a 10% margin of J
  # newabund<- rlnorm(S, meanlog = mm, sdlog = 1)
  # sum(newabund)
  #
  # if (sum(newabund) < J){
  #     repeat{
  #       mm <- mm + .01
  #       newabund<- rlnorm(S, meanlog = mm, sdlog = 1)
  #       if(sum(newabund) < J+(J*.1) && sum(newabund) > J-(J*.1)) break
  #     }
  #   }



  ### if species have the same abundance they are randomly placed in different ranks (within the original rank zone)
  speciesrank=rank(nmatab, ties.method = "random")
  rankab<- newabund[speciesrank]

  ###replacing the matrix list abundances

    for (k in 1:length(matlist[[o]])){
      if (length(rankab)> 0){
        newmat[[o]][[k]] <-rankab[k]
      }
    }

    rankab <- rankab[(-c(1:k))]
  }

  return(newmat)
}


