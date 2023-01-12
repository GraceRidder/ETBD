#!/usr/bin/Rscript

# Sourcing functions ----
source("FC_namesSpecies.R")
source("FC_grow_tree.R")
source("FC_siteSAD.R")
source("FC_Migrate.R")
source("FC_DeleteExtinct.R")
source("FC_DeleteDups.R")
source("FC_Collection_Sep27.R")
source("FC_dist.R")

# Sourcing the simulation ----
#source("ETBDsimOCT_06_2022SPC.R")
source("ETBDsimOCT_24_2022SPC.R")

# #run the simulation multiple times in parallel
trials = 5

repsim <- function(trials){
reslist <- list()
highsp <- data.frame(row.names = 1:12)
highab <-  data.frame(row.names = 1:12)
JmaxV <- c(4000, 8000, 12000, 16000, 20000, 24000, 28000,32000, 36000,40000, 45000,50000)

for (u in 1:trials){
  res1=ETBDspace(t = 1000, pallo = .0, psymp = 0.12, SAD = T, watchgrow = F, SADmarg = .1, siteN = 12, probleave = .0, Jmax = JmaxV, isGrid= F)
  s1 <- c()
  a1 <- c()
  reslist[[u]] <- res1
  for (x in 1:length(JmaxV)){
    s1[x] <- length(unique(unmatrixlist(res1$mig[[length(res1$mig)]][x])))
    a1[x] <- sum(unlist(res1$mig[[length(res1$mig)]][x]))
  }
  highsp <- cbind(highsp, s1)
  highab <- cbind(highab, a1)
}

#result object
results=list(highsp = highsp,
             highab = highab,
             reslist = reslist
)

#result
return(results)

}

rs= mcmapply(repsim, 5, SIMPLIFY = F, mc.cores = 12)

save(rs, file = "rs_5rep_12core_1000t.RData")






#
# lowsp <- data.frame(row.names = 1:20)
# lowab <- data.frame(row.names = 1:20)
# for (u in 1:2){
#   res2=ETBDspace(t = 10, pallo = .0, psymp = 0.08, SAD = T, watchgrow = F, SADmarg = .1, siteN = 20, probleave = .0, Jmax = JmaxV, isGrid= F)
#
#   s2 <- c()
#   a2 <- c()
#   for (x in 1:length(JmaxV)){
#     s2[x] <- length(unique(unmatrixlist(res2$mig[[length(res2$mig)]][x])))
#     a2[x] <- sum(unlist(res2$mig[[length(res2$mig)]][x]))
#   }
#
#   lowsp <- cbind(lowsp, s2)
#   lowab <- cbind(lowab, a2)
#
# }
#
#
#
#
# save(lowsp, file = "lowsp.RData")
# save(lowsp, file = "lowab.RData")
#
#
#
#
# #
# #
# #
# #
