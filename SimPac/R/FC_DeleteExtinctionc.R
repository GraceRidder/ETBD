
##deleting species

DeleteExtinct <- function(mo){

#Saveing exctinct species locations
d <- c()
l <- c()
for( o in 1:length(mo)){
  for (k in 1:length(mo[[o]])){
    if (length(mo[[o]]>0)){
      if(!is.na(mo[[o]][k])){
      if( mo[[o]][k] == 0){
        d <- cbind(d,k)
        l <-cbind(l,o)
      }
    }
  }
  }
}


dl <- rbind(d, l)

dl
#reverse order so I can remove them from the list without changing list order
if (length(l)>1){
  odl<- dl[,order(dl[1,], decreasing = T)]
} else {
  odl <- dl
}


#correct syntax to avoid empty matrix naming issues
if (length(l)>0) {
  for (k in 1:length(l)) {
    if (length(mo[[odl[2, k]]]) > 0) {
      oo <- mo[[odl[2, ][k]]] [-odl[1, ][k],  , drop = FALSE]
      mo[[odl[2, ][k]]] <- as.matrix(oo)
    }
  }
}




if (NA %in% unlist(mo)) {
  message(paste("Problem with deleting the extincts (or an empty matrix)"))
}

return(mo)
}
