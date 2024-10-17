


AlloSpec <- function(matrix_list55, allo_sp, old_sp, Atrip2, splitparm, siteN){

matrix_list4 <- matrix_list55
bud <- T
split <- F


#update the speciating tips
trim <- Atrip2

if (length(unlist(allo_sp)) > 0) {
  g <- 0
  tri <- list()
  triw <- list()
  for (o in 1:length(trim)) {
    tri[[g + 1]] <- trim[[o]][1]
    triw[[g + 1]] <- trim[[o]][2]
    g <- g + 1
  }


  ##for budding speciation one branch has same abundance and new branch has 10% of original
  if (bud) {
    i  <- 1
    fax <- list()
    for (o in 1:length(old_sp)) {
      if (length(old_sp[[o]]) != 0) {
        #update species names in sizes table
        m = match(old_sp[[o]], row.names(matrix_list4[[o]]))
        for (k in 1:length(m)) {
          row.names(matrix_list4[[o]])[m][k] = paste(tri[i], sep = "")
          fax[i] <- (matrix_list4[[o]])[m][k]
          i <- i + 1
        }
      }
    }
  }


  ##for splizting speciation 10% is subtracted from original and new branch is 10% of original
  # if (split) {
  #   i  <- 1
  #   fax <- list()
  #   for (o in 1:length(old_sp)) {
  #     if (length(old_sp[[o]]) != 0) {
  #       #update species names in sizes table
  #       m = match(old_sp[[o]], row.names(matrix_list4[[o]]))
  #       for (k in 1:length(m)) {
  #         row.names(matrix_list4[[o]])[m][k] = paste(tri[i], sep = "")
  #         fax[i] <- (matrix_list4[[o]])[m][k]
  #         faax <- (matrix_list4[[o]])[m][k]
  #         faxtax <- as.numeric(faax) * .1
  #         matrix_list4[[o]][m][k] <-
  #           matrix_list4[[o]][m][k] - faxtax
  #         i <- i + 1
  #       }
  #     }
  #   }
  # }


  #10% of parent population abundance
  flop <- as.matrix(as.numeric(fax) * splitparm)


  ### pop is new species sizes and the new names
  pop <- allo_sp


  i <- 1
  for (o in 1:length(allo_sp)) {
    if (length(allo_sp[[o]]) >= 1) {
      for (k in 1:length(allo_sp[[o]])) {
        pop[[o]][k] <- flop[[i]]
        i <- i + 1
      }
    }
    pop[[o]] <- matrix(as.numeric(pop[[o]]))
  }


  #### altering original names of species
  i <- 1
  for (o in 1:length(allo_sp)) {
    if (length(allo_sp[[o]]) >= 1) {
      m = 1:length(allo_sp[[o]])
      for (k in 1:length(m)) {
        rownames(pop[[o]])[m][k] = paste(triw[i], sep = "")
        i <- i + 1
      }
    }
  }

  #### binding new species onto matrix
  morto <- list()
  for (o in 1:length(siteN)) {
    for (h in 1:length(pop[[o]])) {
      mart <- rbind(matrix_list4[[o]], pop[[o]])
      morto[[o]] <- mart
    }
  }



  # Fix empty row names
  for (o in 1:length(siteN)) {
    for (k in length(morto[[o]]))
      if (length(rownames(morto[[o]])) < 1) {
        rownames(morto[[o]]) <- rownames(matrix_list4[[o]])
      }
  }
  matrix_list5 <- morto
} else {
  matrix_list5 <- matrix_list4
}


if (NA %in% unlist(matrix_list5)) {
  message(
    "Problem with allopatric speciaiton"
  )
  print(matrix_list5)
}

return(matrix_list5)
}




