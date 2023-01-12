
#setwd("~/Desktop/grace/SimCo/")

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
source("ETBDsimOCT_18_2022SPC.R")


res1 = ETBDspace(
  t = 2000,
  pallo = .0,
  psymp = .12,
  SAD = T,
  watchgrow = F,
  SADmarg = .1,
  siteN = 20,
  probleave = .0,
  Jmax = JmaxV,
  isGrid = F
)

save(res1, file = "BigRun2000.RData")








