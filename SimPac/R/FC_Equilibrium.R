
dd <- append(dd, summat)

dd <- sample(1:150, 900, replace = T)

i = 400
ss <- dd[i:(i+30)]

sd(ss)


TS = 50 ## time to reach equilibrium
SDd = 35 ## max deviation allowed for equilibrium

if (length(dd) > TS) {
  for (i in 1:(length(dd)-TS)) {
    ss <- dd[i:(i + TS)]
    sdss <- sd(ss)
    if (sdss < SDd) {
      print(paste("equilibrium reached SD:", sdss, "time:", i))
      break

      end

    }
  }
}


