#Welcome to the migration function  

  Migrate <- function(matrixlist, dist, mig_percent, siteN){
  mote <- matrixlist
  matrix_list1 <- matrixlist
  
  ##samples site locations based on dist(X)
  mig_site <- c()
  mig_loc <- c() 
  for (k in 1:length(mote)) {
    for (o in 1:length(mote[[k]])) {
      if( length(mote[[k]]) > 0){
        probV = dist[k,]
        mote[[k]] <- cbind(sample(1:length(siteN),length(mote[[k]]),replace=TRUE, prob = probV)) 
        mig_site <- append(mig_site, k)
        mig_loc <- append(mig_loc, o)
      }
    }
  } 
  
  ##turns species that don't move to zeros 
  for (o in 1:length(siteN)){
    for (k in 1:length(mote[[o]])){
      if (length(mote[[o]]) > 0){
        if(mote[[o]][k] == o){
          mote[[o]][k] = 0
        }
      }
    }
  }
  

 
  #organizing some information for later 
  new_site <- unlist(mote)
  mig_size <- unlist(matrix_list1)
  mig_sp <- unmatrixlist(matrix_list1)
  hg <- cbind(mig_site, mig_loc, new_site, mig_sp, mig_size)
  hgloc <- subset(hg, new_site != "0")
  
  ## I use hgloc later 
  hgloc <- subset(hgloc, !duplicated(hgloc[,4]))
  new_site <- as.numeric(hgloc[,3])
  mig_size <- as.numeric(hgloc[,5])
  mig_sp <- hgloc[,4]
  mig_site <- as.numeric(hgloc[,1])
  mig_loc <- as.numeric(hgloc[,2])
  hgloc

  ##create a blank matrix list 
  blank1 <-matrix(nrow=0,ncol=1)
  blank <- list()
  for (o in 1:length(siteN)){
    blank[[o]] <- blank1
  }
  
  between <- blank
  if (length(new_site>0)){
    ##create lists w/ species names 
    between <- blank
    letween <- blank
    for (b in 1:length(new_site)){
      if (new_site[[b]] != 0){
        between[[new_site[b]]] <- rbind(between[[new_site[b]]], mig_sp[b]) 
        letween[[new_site[b]]] <- rbind(letween[[new_site[b]]], mig_size[b]) 
      }
    }
    
    # put them together 
    for (o in 1:length(between)){
      rownames(letween[[o]]) <- between[[o]][,1]
    }    
    
    #Make size half of the initial population 
    for ( o in 1:length(siteN)){
      for ( k in 1:length(letween[[o]])){
        if (length(letween[[o]]>0)){
          letween[[o]][[k]] <-letween[[o]][[k]]*mig_percent 
        }
      }
    }
    
    ##bind new population to original matrix 
    more <- list()
    for (o in 1:length(siteN)){
      mart <- rbind(matrix_list1[[o]], letween[[o]])
      more[[o]]<-mart
    }
    
    
    ### to drop duplicates and add size to rowname in individual matricies. 
    temp_size <- list()
    temp_name <- list()
    
    for (o in 1:length(siteN)){
      temp_size[[o]] <- more[[o]][duplicated(row.names(more[[o]]))]
      temp_name[[o]] <- row.names(more[[o]])[duplicated(row.names(more[[o]]))]
      if (length(which(duplicated(row.names(more[[o]])))) > 0 ){
        more[[o]] <- as.matrix(more[[o]][-which(duplicated(row.names(more[[o]]))),])
      }
    }
    
    #fixing bug with small sites 
    for ( o in 1:length(siteN)){
      for ( k in 1:length(more[[o]])){
        if(length(more[[o]]) == 1){
          row.names(more[[o]])[k] <- row.names(more[[o]])[k]
        }
      }
    }
    
    if (length(temp_name) > 0 ){
      for (o in 1:length(temp_name)){
        for (k in 1:length(more[[o]])){
          for(b in 1:length(temp_name[[o]])){
            if (length(temp_name[[o]]) > 0 ){
              if (temp_name[[o]][[b]] == row.names(more[[o]])[k]){
                more[[o]][[k]]  <- more[[o]][[k]] +temp_size[[o]][[b]]
              }
            } 
          }
        }
      } 
    }
    
    
    
    
  } else {
    more <- matrix_list1
  }
  if (NA %in% unlist(more)) {
    warning(paste("Problem with migration"))
  }
  
  results = list(matrixlist = more,
                 data = hgloc, 
                 allo = between)
    return(results)
  }



 
 
 
    