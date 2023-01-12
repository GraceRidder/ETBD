###function which finds new species names 


## x is the speciating species 
## y is the matrix list the speciating species are from 
new.sp.names <- function(speciesMatrix, MatrixList, pine){
  
  
  #temp hold all unique species
  tampa <- unlist(speciesMatrix)
  #tampa it just a list of spec
  temp <- unique(tampa)
  ##counts how many speciating species exist a cross the sites 
  xlist <- list()
  if (length(temp) > 0){
    for (f in 1:length(temp)){
      xlist[f] <- list(f)
      for( o in 1:length(MatrixList)){
        if( temp[f] %in% row.names(MatrixList[[o]])){
          xlist[[f]] <- append(xlist[[f]], temp[f])
        }
      } 
    }
  }  
  
  ##counts the number of speciating for each species 
  tlist <- list()
  if (length(temp)> 0){
    for (f in 1:length(temp)){
      tlist[f] <- list(f)
      for (j in 1:length(tampa)){
        if( tampa[j] == temp[f]){
          tlist[[f]] <- append(tlist[[f]], tampa[j])
        }
      } 
    }
  }  
  
  if (length(temp)>0){
    
    trip2 <- list()
    if( length(temp)> 0){
      for (f in 1:length(temp)){
        #brlen <- (length(tlist[[f]])-1)*2
        b <- 2
        trip <- c()
        trippy <- c()
        s <-paste(temp[f],ten,sep="")    #makes a list of all possible existing species names (only through)
        
        for (i in 1:length(s)){
          if (s[i] %in% pine$tip.label){   #finds the last species and saves the next 2 possible species 
            trip <- c(s[i+1], s[i+2])
            
          }
        }
        
        if (length(trip) < 2){
          trip <- c(s[1], s[2])
        }
        
        jac <- c()
        for (z in 1:length(tlist[[f]])){
          if (temp[f] == tlist[[f]][z]){
            jac <- append(jac, temp[f])
          }
        }
        
        while(length(trip) < length(jac)*2){
          for (i in 1:length(s)){
            if (s[i] %in% trip){  
              trippy <- c(s[i+1], s[i+2])
            }
          }
          trip <- append(trip, trippy)
          b <- b+2
        }
        
        if (length(xlist[[f]]) > length(tlist[[f]])){
          trip <- append(trip, paste(temp[f], sep=""))
          b <- b + 1
        }
        
        trip2[[f]] <-trip
        
      }
    }
  }
  return(trip2)
}


