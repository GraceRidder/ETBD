

makeLineDomain <- function(siteN, probleave){

probstay = 1-probleave
x <-matrix(data = rep(0, siteN*siteN), nrow = siteN, byrow = TRUE)

for (i in 1:length(x[,1])){
  x[i,i] <- probstay
  if (i>1 && i<length(x[1,])){
  x[i,i-1] <- probleave/2
  x[i,i+1] <- probleave/2
  }
}

x[1,2] <- probleave
x[length(x[,1]),length(x[,1])-1] <-probleave


return(x)

}


## make grid domain
makeGridDomain <- function(siteN, probleave){
probstay = 1-probleave
x <- matrix(1:(siteN*siteN), nrow = siteN, byrow = T)

if (siteN %% sqrt(siteN) != 0) {
  warning(paste("site number not a square"))
}

####   [] [] []          [] [] [] []
####   [] [] []          [] [] [] []
####   [] [] []          [] [] [] []
###                      [] [] [] []
###


rw = sqrt(siteN)
x <- matrix(rep(0, (siteN*siteN)), nrow = siteN, byrow = T)

for (i in 1:length(x[,1])){
  x[i,i] <- probstay

  ##first row
  if (i>1 && i<rw){
    x[i,i-1] <- probleave/3
    x[i,i+1] <- probleave/3
    }
  ##first row  middles
  if(i<rw){
    x[i,i+rw] <- probleave/3
  }

  ##last row
  if (i<length(x[,1]) && i >length(x[,1])-rw){
    x[i,i-1] <- probleave/3
    x[i,i+1] <- probleave/3
  }
  ##last row middles
    if(i> length(x[1,])-rw){
      x[i,i-rw] <- probleave/3
    }
    ##top right
    if(i==rw){
     x[i,i-1] <- probleave/2
     x[i,i+rw] <- probleave/2
    }
     ##bottom left
     if ( i==(length(x[,1])-(rw-1)) ){
       x[i,i-1] <- 0
       x[i,i+1] <- probleave/2
       x[i,i-rw] <- probleave/2
     }
   ##right edge
  right <- c()
  for(l in 1:(siteN/rw)){
    right<- append(right, rw*l)
  }

  right <- right[-c(1)]
  right <-  right[-c(length(right))]

  if (i %in% right){
    x[i,i-1] <- probleave/3
    x[i,i-rw] <- probleave/3
    x[i,i+rw] <- probleave/3
  }

  ##left edge
  left <- c()
  for(l in 1:(siteN/rw)){
    left<- append(left, (rw*l)+1)
  }

  left <-  left[-c(length(left))]
  left <- left[-c(length(left))]

  if (i %in% left){
    x[i,i+1] <- probleave/3
    x[i,i-rw] <- probleave/3
    x[i,i+rw] <- probleave/3
  }
    ##middle squares
  for ( k in 1:length(right)) {
    if (  i > left[k] && i < right[k]){
      x[i,i-1] <- probleave/4
      x[i,i+1] <- probleave/4
      x[i,i-rw] <- probleave/4
      x[i,i+rw] <- probleave/4
    }
  }
    }


x[1,2] <- probleave/2
x[1,rw+1] <- probleave/2

x[siteN,siteN-1] <- probleave/2
x[siteN,siteN-(rw)] <- probleave/2


return(x)
}

