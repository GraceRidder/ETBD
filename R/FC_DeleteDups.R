#deleting duplicated allopatrically speciating species

DeleteDups <- function(matrixlist){
  
cafe <- matrixlist 
i <- 1
dup <- duplicated(unlist(matrixlist))
for (o in 1:length(matrixlist)){
  for (k in 1:length(matrixlist[[o]])){
    if (length(matrixlist[[o]]>0)){
      cafe[[o]][k] <- dup[i]
      i <- i+1
    }
  }
}

ll <- c()
dd <- c()
for (o in 1:length(matrixlist)){
  for (k in 1:length(matrixlist[[o]])){
    if (length(matrixlist[[o]]>0)){
      if (cafe[[o]][k]== TRUE){
        dd <- cbind(dd,k)
        ll <-cbind(ll,o)
      }
    }
  }
}

ddll <- rbind(dd, ll)

#reverse order so I can remove them from the list without changing list order
if (length(ll)>1){
  oddll<- ddll[,order(ddll[1,], decreasing = T)]
} else {
  oddll <- ddll
}


#correct syntax to avoid empty matrix naming issue 
if (length(ll)>0){
  for( k in 1:length(ll)){
    oo <- matrixlist[[oddll[2,][k]]] [-oddll[1,][k],  ,drop = FALSE]
    matrixlist[[oddll[2,][k]]] <- as.matrix(oo)
  }
}

return(matrixlist)
}







