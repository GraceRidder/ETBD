# Setting up workspace ----

# If you are using Base R, this sets the file path to wherever this script is:
#setwd(getSrcDirectory()[1])

# On R-Studio, use:
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Sourcing functions ----
# source("FC_namesSpecies.R")
# source("FC_grow_tree.R")
# source("FC_SAD.R")
# source("ETBD_Initialization_Jun13.R")
# source("FC_Migrate.R")
# source("FC_DeleteExtinct.R")
# source("FC_DeleteDups.R")





ETBDspace = function(initialtree,
                     #initial phylogeny
                     t = 10,
                     #no timesteps
                     # EXCALC = 1/(0.37^1.1)*exp(-1.1*matrix_list5[[o]][,1]),
                     ### EXCALC = exp(-.2*matrix_list5[[o]][,1]),
                     #EXCALC = exp(-.1*matrix_list5[[o]][,1]),
                    # EXCALC = .0,
                     E0 = 1.1,
                     #extinction convexity
                     Ec = 1 / (0.37 ^ 1.1),
                     #extinction linear scaling ext probability=Ec*exp(-E0*N)
                     Jmax = 100000 ,
                     #uper bound of individuals
                     JmaxV =  c(
                       1500,
                       2000,
                       2500,
                       3000,
                       3500,
                       4000,
                       8000,
                       12000,
                       16000,
                       20000,
                       24000,
                       28000,
                       32000,
                       36000,
                       40000,
                       45000,
                       50000,
                       55000,
                       60000,
                       65000
                     ),
                     MIG = T,
                     split = T,
                     ## speciation type can be either "split" or "bud"
                     bud = F,
                     siteN = 20,
                     probleave = .1,
                     SAD = T,
                     #psymp = ((matrix_list1[[o]][,1])/JmaxV[o])^.4,
                     psymp = JmaxV[o] / 100000,
                     pallo = .3,
                     watchgrow = T,
                     mig_percent = .3,
                     SADmarg = .1,
                     isGrid = F)


  {{
    #monitor objects
    extincttotal = NULL
    migration = NULL
    trip2 <- NULL
    trees = list()
    mig = list()
    extinctsp = list()
    allop = list()
    symp = list()
    allotrip = list()
    symptrip = list()
    exsp = list()

    ##run these to run a time step individually
    # ipa = 1
    # psymp = .5
    # pallo = .5
    #  MIG = T
    #  split = T
    #  bud = F
    # allopatric = T
    # siteN = 20
    # probleave = .3
    # SAD = T
    # watchgrow = T
    # isGrid = F
    # mig_percent = .3
    # SADmarg = .1


    #Define the spatial domain
    if (isGrid) {
      dist <- makeGridDomain(siteN, probleave)
    } else{
      dist <- makeLineDomain(siteN, probleave)
    }


    #initial tree with 100 dist
    {
      #initial species sizes
      siteN = 1:siteN
      grow <- c(4, 4, 4, 4, 4)
      initialsize = 100
      names = as.factor(c(paste(
        "t", str_pad(siteN, 3, pad = "0"), sep = ""
      )))
      names2 = as.factor(c(paste(
        "t", str_pad(siteN, 3, pad = "0"), sep =
          ""
      )))
      tree = as.phylo(~ names)
      tree$edge.length = rep(1, length(siteN))

      ITtree <- tree

      matrix_list <- list()
      for (s in 1:length(siteN)) {
        q <- matrix(1)
        row.names(q) = tree$tip.label[s]
        matrix_list[[s]] = q
      }

      ma <- list()
      for (s in 1:length(siteN)) {
        q <- matrix(
          c(1),
          nrow = 1,
          ncol = 1,
          byrow = FALSE,
          dimnames = list(c("start0"),
                          c("sizes"))
        )
        ma[[s]] = q
      }

      matrix_list0 <- matrix_list
      test1 <- ma

      for (o in 1:length(matrix_list0)) {
        matrix_list0[[o]] <- matrix_list0[[o]] + initialsize - 1
      }
    }

    temptree <- ITtree
    matrix_temp <- matrix_list0
    matrix_list6 <- matrix_list0
    blank1 <- matrix(nrow = 0, ncol = 1)
    blank <- list()
    for (o in 1:length(siteN)) {
      blank[[o]] <- blank1
    }

     #matrix_list5 <-  matrix_list0
    # EXCALC = 1/(0.37^1.1)*exp(-1.1*matrix_list5[[o]][,1])
    #psymp = JmaxV[o]/40000



    for (ipa in 1:t)
      {
        #start
        matrix_list0 <- matrix_list6
        matrix_list0
        ##deleting extinct species from matrix list 6 from previous step
        matrix_list1 <- DeleteExtinct(matrix_list0, siteN)

        if (MIG) {
          ###Migration function: input a matrix list, a probability matrix, and the percent of population that leaves
          migrated <- Migrate(matrix_list1, dist, mig_percent, siteN)

          ###matrix list with species migrated but not speciated
          matrix_list2 <- migrated$matrixlist

          ###Data we need later on where the species came from, where they went, names, etc
          hgloc <- migrated$data
          allo <- migrated$allo

        } else {
          matrix_list2 <- matrix_list1
        }

        ## counting ten time steps
        is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
        B <- ipa/10
        if (is.wholenumber(B)){
          print(paste('time continues on...', ipa))
        }


        matrix_list3 <- matrix_list2

        #sympatrically speciating species



          {
            speciatinglog = list()
            for (o in 1:length(siteN)) {
              spec = as.logical(rbinom(length(matrix_list1[[o]]), 1, psymp))  ##probability of sympatric speciation psymp
              speciatinglog[[o]] = spec
            }


            speciating = list()
            for (o in 1:length(siteN)) {
              spe = setdiff(rownames(matrix_list1[[o]])[speciatinglog[[o]]], extincttotal)
              speciating[[o]] = spe
            }
            symp_sp <- list()
            for (o in 1:length(speciating)) {
              speciatin <- subset(speciating[[o]], speciating[[o]] != "1")
              speciatin <- as.matrix(speciatin)
              symp_sp[[o]] <- speciatin
            }
          }




          symp_sp <- DeleteDups(symp_sp)

          #allopatrically speciating species
          {
            kelplog = list()
            for (o in 1:length(siteN)) {
              kelp = as.logical(rbinom(length(allo[[o]]), 1, pallo))  ## allopatrically speciaiting currently palla (species that migrate have palla probability they speciaite)
              kelplog[[o]] = kelp
            }
            kelping = list()
            for (o in 1:length(siteN)) {
              spe = setdiff(allo[[o]][kelplog[[o]]], extincttotal)
              kelping[[o]] = spe
            }

            allo_sp <- list()
            for (o in 1:length(kelping)) {
              krep <- subset(kelping[[o]], kelping[[o]] != "1")
              krep <- as.matrix(krep)
              allo_sp[[o]] <- krep
            }
          }
          allo_sp <- DeleteDups(allo_sp)

          #making sure allo and symp don't overlap i.e. species cannot allopatrically and sympatrically speciate even if they are different populations
          unsymp <- unlist(symp_sp)
          {
            ll <- c()
            dd <- c()
            for (o in 1:length(siteN)) {
              for (k in 1:length(allo_sp[[o]])) {
                if (length(allo_sp[[o]]) != 0) {
                  if (length(unsymp) > 0) {
                    if (allo_sp[[o]][k] %in% unsymp) {
                      dd <- cbind(dd, k)
                      ll <- cbind(ll, o)
                    }
                  }
                }
              }
            }

            ddll <- rbind(dd, ll)

            #reverse order so I can remove them from the list without changing list order
            if (length(ll) > 1) {
              oddll <- ddll[, order(ddll[1, ], decreasing = T)]
            } else {
              oddll <- ddll
            }

            #correct syntax to avoid empty matrix naming issues
            if (length(ll) > 0) {
              for (k in 1:length(ll)) {
                oo <- allo_sp[[oddll[2, ][k]]] [-oddll[1, ][k],  , drop = FALSE]
                allo_sp[[oddll[2, ][k]]] <- as.matrix(oo)
              }
            }



            if (NA %in% unlist(matrix_list3)) {
              warning(paste("NA in matrixlist3", ipa))
            }

          }

          #finding the species names for the allos
          temptree <- tree
          allo <- allo_sp
          trip2 <- new.sp.names(allo_sp, matrix_list1, temptree)

          #saving allotrip for later
          if (length(unlist(allo_sp)) > 0) {
            allotrip <- trip2
          } else {
            allotrip <- list()
          }


          #update the allopatric names in matrix
          temp <- unique(unlist(allo_sp))

          if (length(temp) > 0) {
            trim <- trip2

            g <- 0
            tri <- list()
            triw <- list()
            for (o in 1:length(trim)) {
              tri[[g + 1]] <- trim[[o]][1]
              triw[[g + 1]] <- trim[[o]][2]
              g <- g + 1
            }

            matrix_list4 <- matrix_list3

            i <- 1
            for (o in 1:length(allo_sp)) {
              if (length(allo_sp[[o]]) > 0) {
                for (k in 1:length(matrix_list3[[o]])) {
                  if (row.names(matrix_list3[[o]])[k] %in% allo_sp[[o]]) {
                    row.names(matrix_list4[[o]])[k] <-  paste(triw[i], sep = "")
                    i <- i + 1
                  }
                }
              }
            }

            ###matrix_list3_4 is now with one new named species and the old species untouched
            matrix_list3_4 <- matrix_list4

            u <- 1
            for (o in 1:length(allo_sp)) {
              for (c in 1:length(hgloc[, 4])) {
                if (hgloc[c, 4] %in% allo_sp[[o]]) {
                  row.names(matrix_list4[[as.numeric(hgloc[c, 1])]])[as.numeric(hgloc[c, 2])] <-
                    tri[u]
                  u <- u + 1
                }
              }
            }

          } else {
            matrix_list4 <- matrix_list3
          }


          if (NA %in% unlist(matrix_list4)) {
            warning(paste("NA in matrixlist4", ipa))
          }


          #######end of allopatric speciation##########

          #growing the tree with allo species
          allotree <-
            grow.tree(allo_sp, matrix_list3_4, temptree, temptree)$tree

          ###add the symp species onto the allotree
          full.tree <-
            grow.tree(symp_sp, matrix_list2, allotree, temptree)$tree

          ###finding names for sympatrically speciating species with homemade function
          trip2 <- NULL
          trip2 <- new.sp.names(symp_sp, matrix_list1, temptree)
          trip2 <- grow.tree(symp_sp, matrix_list2, allotree, temptree)$trip2

          #temp hold all unique species
          temp <- unique(unlist(symp_sp))
          symptrip <- trip2

          #update the speciating tips
          trim <- trip2

          if (length(temp) > 0) {
            g <- 0
            tri <- list()
            triw <- list()
            for (o in 1:length(trim)) {
              tri[[g + 1]] <- trim[[o]][1]
              triw[[g + 1]] <- trim[[o]][2]
              g <- g + 1
            }


            ##for budding speciation one branch has same abundance and new branch has 10% of origianl
            if (bud) {
              i  <- 1
              fax <- list()
              for (o in 1:length(symp_sp)) {
                if (length(symp_sp[[o]]) != 0) {
                  #update species names in sizes table
                  m = match(symp_sp[[o]], row.names(matrix_list4[[o]]))
                  for (k in 1:length(m)) {
                    row.names(matrix_list4[[o]])[m][k] = paste(tri[i], sep = "")
                    fax[i] <- (matrix_list4[[o]])[m][k]
                    i <- i + 1
                  }
                }
              }
            }

            ##for splitting speciation 10% is subtracted from original and new branch is 10% of original
            if (split) {
              i  <- 1
              fax <- list()
              for (o in 1:length(symp_sp)) {
                if (length(symp_sp[[o]]) != 0) {
                  #update species names in sizes table
                  m = match(symp_sp[[o]], row.names(matrix_list4[[o]]))
                  for (k in 1:length(m)) {
                    row.names(matrix_list4[[o]])[m][k] = paste(tri[i], sep = "")
                    fax[i] <- (matrix_list4[[o]])[m][k]
                    faax <- (matrix_list4[[o]])[m][k]
                    faxtax <- as.numeric(faax) * 0.1
                    matrix_list4[[o]][m][k] <-
                      matrix_list4[[o]][m][k] - faxtax
                    i <- i + 1
                  }
                }
              }
            }


            #10% of parent population abundance
            flop <- as.matrix(as.numeric(fax) * 0.1)


            ### pop is new species sizes and the new names
            pop <- symp_sp


            i <- 1
            for (o in 1:length(symp_sp)) {
              if (length(symp_sp[[o]]) >= 1) {
                for (k in 1:length(symp_sp[[o]])) {
                  pop[[o]][k] <- flop[[i]]
                  i <- i + 1
                }
              }
              pop[[o]] <- matrix(as.numeric(pop[[o]]))
            }


            i <- 1
            for (o in 1:length(symp_sp)) {
              if (length(symp_sp[[o]]) >= 1) {
                m = 1:length(symp_sp[[o]])
                for (k in 1:length(m)) {
                  rownames(pop[[o]])[m][k] = paste(triw[i], sep = "")
                  i <- i + 1
                }
              }
            }


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
            warning(
              paste(
                "NA in matrixlist5: problem with adding sympatric species to matrix list",
                ipa
              )
            )
          }



          ########updating survivors###########

          ##list of speciating species
          specspec <- append(unlist(allo_sp), unlist(symp_sp))

          ##list of new species that have already grown in tree
          grownspec <- unlist(unique(append(symptrip, allotrip)))
          allspec <-  unmatrixlist(matrix_list5)

          ##finding species that didn't speciate or grow or are extinct.
          sop = setdiff(allspec, specspec)
          sop2 = setdiff(sop, grownspec)
          sop3 = setdiff(sop2, extincttotal)
          ma <- test1
          ma[[1]] <- unique(sop3)

          #Updating survivors in tree
          tree <- full.tree
          if (length(sop3 > 1)) {
            for (o in 1:length(siteN)) {
              for (k in 1:length(ma[[o]])) {
                if (ma[[o]][k] != 1) {
                  tree = surviveatx(tree, ma[[o]][k])
                }
              }
            }
          }

          if (watchgrow) {
            plot(tree, cex = .5)
          }



          ### RANKS ABUNDANCES AND DRAWS FROM SAD Fishers log series distribution
          if (SAD) {
            if (length(unmatrixlist(matrix_list5)) > 5) {
              xx <- MakeSAD(matrix_list5, SADmarg, JmaxV)
              matrix_list5 <- xx
            }
          }

          if (NA %in% unlist(matrix_list5)) {
            warning(paste(
              "NA in matrixlist5: problem with the SAD rank setting",
              ipa
            ))
          }

          ####### extinction #########
          {


            #calculate extinction probability
            etip = list()
            for (o in 1:length(siteN)) {
              extinctionp = 1 / (0.37 ^ 1.1) * exp(-1.1 * matrix_list5[[o]][, 1])
              etip[[o]] <- extinctionp
            }



            for (o in 1:length(siteN)) {
              for (i in 1:length(etip[[o]])) {
                if (etip[[o]][i] > 1) {
                  etip[[o]][i] <- 1
                }
              }
            }


            #as logical...
            etipp = list()
            for (o in 1:length(siteN)) {
              extlog = as.logical(rbinom(length(matrix_list5[[o]][, 1]), 1, etip[[o]]))
              etipp[[o]] <- extlog
            }


            #name and position of extinct
            extinct = list()
            for (o in 1:length(siteN)) {
              if (length(matrix_list5[[o]] > 0)) {
                ex <- row.names(matrix_list5[[o]])[etipp[[o]]]
                extinct[[o]] <- as.matrix(ex)
              }
            }

            extinctx <- extinct

            ext <- c()
            for (o in 1:length(extinct)) {
              if (length(extinct[[o]]) > 0) {
                ext <- append(ext, (extinct[[o]]))
              }
            }

            matrix_list6 <- matrix_list5

            if (NA %in% unlist(matrix_list6)) {
              warning(paste(
                "NA in matrixlist6: problem with the extinction",
                ipa
              ))
            }


            for (o in 1:length(extinctx)) {
              if (length(extinctx) > 0) {
                if (length(extinctx[[o]]) > 0) {
                  #prune the extinct from table
                  m = match(extinctx[[o]], row.names(matrix_list5[[o]]))
                  matrix_list6[[o]][m, 1] = 0
                }
              }
            }


            extincttotal = c(extincttotal, ext)



            summat <- sum(unlist(matrix_list6))

            if (sum(summat) == 0){
              stop('everything is dead')
              }



# break if a site goes extinct
# break if we reach equilibrium
          #300 time steps within some range 10 species?







          #   ##growing things based on grow
          #   for (o in 1:length(siteN)) {
          #     if (length(matrix_list6[[o]]) > 0) {
          #       for (k in 1:length(matrix_list6[[o]])) {
          #         matrix_list6[[o]][[k]] <-
          #           matrix_list6[[o]][[k]] + matrix_list6[[o]][[k]] / grow[o]
          #       }
          #     }
          #   }
          # }


          # ### RANKS ABUNDANCES AND DRAWS FROM SAD Fishers log series distribution
          # if (SAD) {
          #   if (length(unmatrixlist(matrix_list5)) > 5) {
          #     xx <- MakeSAD(matrix_list6, SADmarg, JmaxV)
          #     matrix_list6 <- xx
          #   }
          # }









          ###### End of simulation
          ####Starts again with matrix_list6 at the start


          #
          #
          # ####everything beyond this point is just to help me catch issues
          # x <- unmatrixlist(matrix_list6)
          #
          # length(unique(x))
          # length(tree$tip.label)
          #
          #
          # ### species in matrix but not in the tree
          # for (o in 1:length(x)){
          #   if (x[o] %!in% tree$tip.label){
          #     print(paste(x[o], "not in tree"))
          #
          #   }
          # }
          #   hecklist <- c()
          #   i <- 1
          # ##species in the tree but not in the martix .............
          # for (o in 1:length(tree$tip.label)){
          #   if (tree$tip.label[o] %!in% x && tree$tip.label[[o]] %!in% extincttotal){
          #     print(paste(tree$tip.label[o], "not in matrix"))
          #     hecklist[[i]] <- append(hecklist, tree$tip.label[o])
          #     i <- 1+i
          #   }
          # }
          #
          #   if(length(hecklist) > 0){
          #   for ( o in 1:length(hecklist)){
          #     if ( hecklist[o] %in% unlist(symptrip)){
          #       print ("in sympatric speciating list")
          #     }
          #   }
          #
          #
          #       for ( o in 1:length(hecklist)){
          #         if ( hecklist[o] %!in% unlist(tri)){
          #           print ("not in tri")
          #         }
          #       }
          #
          #
          #       for ( o in 1:length(hecklist)){
          #         if ( hecklist[o] %in% unmatrixlist(pop)){
          #           print ("in pop")
          #         }
          #       }
          #
          #   for ( o in 1:length(hecklist)){
          #     if ( hecklist[o] %in% unlist(allotrip)){
          #       print ("in allopatric speciating list")
          #     }
          #   }
          # }
          #
          #     ##species in the tree but not in the martix .............
          #     for (o in 1:length(tree$tip.label)){
          #       if (tree$tip.label[o] %!in% x && tree$tip.label[[o]] %!in% extincttotal){
          #         print(paste(tree$tip.label[o], "not in matrix"))
          #       }
          #     }
          #
          #
          #


        }

        #monitors of sizes and trees
        extinctsp[[ipa]] = ext
        trees[[ipa]] = tree
        mig[[ipa]] = matrix_list6
        allop[[ipa]] = allo_sp
        symp[[ipa]] = symp_sp
        allotrip[[ipa]] = allotrip
        symptrip[[ipa]] = symptrip
        exsp[[ipa]] = extinct
      }

      return(
        list(
          tree = tree,
          #final tree
          trees = trees,
          #all trees by timeslice
          matrix_list = matrix_list6,
          mig = mig,
          extincts = extinctsp,
          allop = allop,
          symp = symp,
          allotrip = allotrip,
          symptrip = symptrip,
          exsp  = exsp

        )
      )

    }
    }


