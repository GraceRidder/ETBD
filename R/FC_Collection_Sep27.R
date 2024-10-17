
library(ape)
library(phytools)
library(geiger)
library(stringr)
library(sads)
library(MASS)


#helper function for solving the geometric SAD
######################################################
geometricBrentsolver=function(S, #number of species
                              J, #total nuber of individuals
                              Nmin=1, #abundance of least abundant species
                              bounds=c(1.000001,100), #bounds of search interval for r
                              ... )
{
  #function to be solved
  fn=function(r) Nmin*((1-r^S)/(1-r))-J

  #looking for more appropriate bounds
  if (fn(bounds[1])>0){
    repeat{
      bounds[1]=1+((bounds[1]-1)/10)
      if(fn(bounds[1])<0) break
    }
  }

  #iterative solver using
  root=uniroot(fn, bounds)
  #result object
  results=list(r=root$root, #r
               Nvect=Nmin*root$root^((S:1)-1), #vector of species abundances
               predict_J=(Nmin*(1-root$root^S)/(1-root$root)), #J prediction based on parameters (to check numerical errors)
               bounds=bounds #final bounds
  )

  #result
  return(results)

}


#helper functions simulating speciating and survivor tips
#with names successively filled by 1 and 2
speciateatx=function(tr, tip){
  #creates a subtree of size 2 and length 1 to be binded into phylogeny
  subtreetobind=ape::rtree(2,br=c(1,1), tip.label = c(paste(tip,2,sep=""),paste(tip,1,sep="")))
  #binding
  tr=ape::bind.tree(tr,subtreetobind,where=which(tr$tip.label==tip))
  return(tr)
}

surviveatx=function(tr, tip){
  #creates a subtree of size 1 and length 1 to be binded to phylogeny
  subtreetobind=ape::rtree(2,br=c(1,1), tip.label = c(paste(tip),paste(tip)))
  onetiptobind=ape::drop.tip(subtreetobind,1)
  #binding
  tr=ape::bind.tree(tr,onetiptobind,where=which(tr$tip.label==tip))
  return(tr)
}


speciateatXX=function(tr, tip){
  #creates a subtree of size 2 and length 1 to be binded into phylogeny
  subtreetobind=ape::rtree(2,br=c(1,1), tip.label = c(paste(tip,sep=""),paste(twix, sep="")))
  #binding
  tr=ape::bind.tree(tr,subtreetobind,where=which(tr$tip.label==tip))
  return(tr)
}


##make a list of possible existing species +1
ten <- c(1:9)

repeat{
  extra <- c()
  for (r in 1:length(ten)){
    extra <- append(extra, paste('0', ten[r], sep = ""))
  }

  ten <- c(1:9, extra)
  if( length(ten) > 130){
    break
  }
}


ten

'%!in%' <- function(x,y)!('%in%'(x,y))
################################



substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

allDuplicated <- function(vec){
  front <- duplicated(vec)
  back <- duplicated(vec, fromLast = TRUE)
  all_dup <- front + back > 0
  return(all_dup)
}


unmatrixlist <- function(x){
  choo <- c()
  for (o in 1:length(x)){
    for (k in 1:length(x[[o]])){
      if (length(x[[o]]) > 0){
        choo <- append(choo, row.names(x[[o]])[k])

      }

    }
  }
  return(choo)
}





trip2 = NULL
