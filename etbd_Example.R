#### Example script on how to run the ETBD simulator


install.packages("~/Desktop/ETBD.tar.xz", repos = NULL)
#### the following will run an example of population size dependent speciation.

res1 = simETBD(
  t =200,
  DIST = "SRS",     ### NO, GEO, SRS, NORM
  JmaxV = 2000,
  ### negative exponential extinction ###
  NegExpEx = F,   ### turns it on
  exparm = -.7,   ### curveature of extinction function
  psymp = .35,    ### constant speciation
  ### negative exponential speciation ###
  ExpSp = T,      ### turns it on
  ExpSpParm = 1.6,  #### strength of the speciation
  constantEX = .005, #### constant extinction
  splitparm = .5, ### speciation symmetry
)

##read the newik
myTree <- ape::read.tree(text = res1$tree)
###how many nodes are on the phylogeny
myTree$Nnode

###plot the phylogeny
plot(myTree, show.tip.label = F)

###plot the phylogeny without extinction
plot(paleotree::dropExtinct(myTree), show.tip.label = F)


###plot the species richness though time
rich <- c()
for (i in 1:length(res1$matrix_lists)){
  rich <- append(rich,length(res1$matrix_lists[i][[1]][[1]][,1]))
}

plot(rich, typ = "l", xlab = "time" ,ylab = "species richness" )

###plot the SAD
hist(res1$matrix_lists[[20]][[1]][,1], main = "SAD", xlab = "abundance")




#### the following will run an example of population size dependent extinction

## WARNING this is much slower because in this case extinction speeds up (and so does speciairton)
## generating many more species and many more general events

res1 = simETBD(
  t =50,
  DIST = "SRS",     ### NO, GEO, SRS, NORM
  JmaxV = 2000,
  ### negative exponential extinction ###
  NegExpEx = T,   ### turns it on
  exparm = -.7,   ### curveature of extinction function
  psymp = .35,    ### constant speciation
  ### negative exponential speciation ###
  ExpSp = F,      ### turns it on
  ExpSpParm = 1.6,  #### strength of the speciation
  constantEX = .005, #### constant extinction
  splitparm = .5, ### speciation symmetry
)

##read the newik
myTree <- ape::read.tree(text = res1$tree)

###how many nodes are on the phylogeny
myTree$Nnode

###plot the phylogeny
plot(myTree, show.tip.label = F)

###plot the phylogeny without extinction
plot(paleotree::dropExtinct(myTree), show.tip.label = F)

###how many nodes are on the phylogeny (without the extinct species)
paleotree::dropExtinct(myTree)$Nnode

###plot the species richness though time
rich <- c()
for (i in 1:length(res1$matrix_lists)) {
  rich <- append(rich, length(res1$matrix_lists[i][[1]][[1]][, 1]))
}
plot(rich,
     typ = "l",
     xlab = "time" ,
     ylab = "species richness")

###plot the SAD
hist(res1$matrix_lists[[20]][[1]][, 1], main = "SAD", xlab = "abundance")


