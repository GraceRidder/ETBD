
totalAB = 0
hablist <- c()
hsplist <- c()
for (x in 1:length(res$matrix_list)){
  totalAB = totalAB + sum(length(res$mig[[1000]][[x]]))
  hablist <- append(hablist, sum(unlist(res$mig[[1000]][[x]])))
  hsplist <- append(hsplist, length(res$mig[[1000]][[x]]))
}

hablist
hsplist

res <- res1
totalAB = 0
hablist <- c()
hsplist <- c()
for (x in 1:length(res$matrix_list)){
  totalAB = totalAB + sum(res$matrix_list[[x]])
  hablist <- append(hablist, sum(res$matrix_list[[x]]))
  hsplist <- append(hsplist, length(res$matrix_list[[x]]))
}

y <- rs[[1]]$highsp
plot(rowMeans(y), typ = "l", col = "red")

x <- rs[[1]]$highab
plot(rowMeans(x), typ = "l", col = "red")

plot(rowMeans(x), rowMeans(y), typ = "l")

#total abundance
totalAB
#Abundance by site
barplot(ablist, main = "Abundance per site", ylab = "Abundance", xlab = "Site (Jmax gradient)")
#species richness per site
plot(hsplist[5:15], typ = "l", col = "red", ylim = c(0,400), main = "Species richness per site", ylab = "Species richness", xlab = "Site (Jmax gradient")
lines(log(lsplist[5:15]), typ = "l", col = "blue")

#species abundance per site
plot(log(hablist[1:20]), typ = "l", col = "red", ylim = c(0,17), main = "Abundance per site", ylab = "Abundance", xlab = "Site (Jmax gradient")
lines(log(lablist[5:15]), typ = "l", col = "blue")

#species richness per site
plot(hablist[5:15], hsplist[5:15], typ = "l", col = "red", ylim = c(0,400), main = "Species and Abundance", xlab = "Abundance", ylab = "Species richness")
lines(lablist[5:15],lsplist[5:15], typ = "l", col = "blue")

#species richness per site
plot(log(hablist[5:15]), log(hsplist[5:15]), typ = "l", col = "red", ylim = c(4,6), main = "Species and Abundance ", xlab = "log(Abundance)", ylab = "log(Species richness)")
lines(log(lablist[5:15]),log(lsplist[5:15]), typ = "l", col = "blue")

res<- res1



datapac <- c()
res <- rs[[1]]$reslist[[5]]


tpac <- c()
apac<- c()
for ( i in 1:length(res$mig)){
  tpac <- append(tpac, length(unlist(res$mig[i][[1]][12])))
  apac <- append(apac, sum(unlist(res$mig[i][[1]][12])))
}



datapac <- cbind(datapac, apac)
ncol(datapac)


aa <- rowMeans(datapac)

plot(aa, typ= "l", main = "Species through time in site 12", xlab = "Time-step", ylab = "Species")

plot(tpac, typ= "l", main = "Species through time in site 12", xlab = "Time-step", ylab = "Species")
plot(apac, typ= "l", main = "Abundance through time in site 11", xlab = "Time-step", ylab = "Abundance")






res<- res1
j <- c()
s <- c()
f = 3
for (x in 1:length(res$mig)){
  j[x] <- sum(unlist(res$mig[[x]][f]))
  s[x] <- length(unlist(res$mig[[x]][f]))
}

dev.off()
plot(s, j, col = "red")
Jmax = 4000
curve(Jmax*(x/(100+x)), from=0, to=max(s), add = T)
Jmax = 12000
curve(Jmax*(x/(100+x)), from=0, to=max(s), add = T)
Jmax = 65000
curve(Jmax*(x/(100+x)), from=0, to=max(s), add = T)
Jmax = 80000
curve(Jmax*(x/(100+x)), from=0, to=max(s), add = T)
Jmax = 100000
curve(Jmax*(x/(100+x)), from=0, to=max(s), add = T)
