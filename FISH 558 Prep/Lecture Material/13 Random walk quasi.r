#Random walk model
#specify number of years, starting population size, probability b of producing new offspring, 
#probability d of dying (and hence probability 1-b-d of neither dying nor producing new offspring)
#assume death implies no reproduction. 
#!!!!Added in a quasiextinction parameter, if N<quasi then N=0 extinction
#Code written by Trevor A. Branch, 21-23 April 2012, tbranch@uw.edu


random.walk <- function(N, b, d, nyears=100, quasiextinction=10) {
  N.vector <- vector(length=nyears)
  if (N < quasiextinction) {   #no point in going through all the years
    N.vector[]<-0  #set all N's to zero
  }
  else {
    N.vector[1] <- N  #first year
    for (yr in 2:nyears) {
      probs <- runif(n=N.vector[yr-1])  #vector of probabilities between 0 and 1
      births <- sum(probs<b)    #probs<b returns vector of T or F, but T=1, and F=0, therefore sum() is the number of T
      deaths <- sum(probs>=b & (probs < b+d))  #if between b and b+d then death
      N.vector[yr] <- N.vector[yr-1]+births-deaths
      if (N.vector[yr] < quasiextinction) {
        N.vector[yr]<-0
      }
    }                
  }
  invisible(N.vector)
}
temp <- cbind(random.walk(N=30,b=0.2,d=0.2,nyears=100),random.walk(N=30,b=0.2,d=0.2,nyears=100),random.walk(N=30,b=0.2,d=0.2,nyears=100),
     random.walk(N=30,b=0.2,d=0.2,nyears=100),random.walk(N=30,b=0.2,d=0.2,nyears=100))
write.csv(temp,file="Tables\\SampleRandomWalksQuasi.csv")

extinction.prob <- function(N.values=seq(1,100,1), b=0.2, d=0.2, nyears=100, nsims=100, quasiextinction=10) {
  p.extinct <- vector(length=length(N.values))
  for (j in 1:length(N.values)) {
    nextinct <- 0
    for (k in 1:nsims) {
      if (random.walk(N=N.values[j],b=b,d=d,nyears=nyears, quasiextinction=quasiextinction)[nyears]==0) {
        nextinct <- nextinct+1
      }
    }
    p.extinct[j] <- nextinct/nsims
  }
  plot(x=N.values,y=p.extinct)
  invisible(cbind(N.values,p.extinct))
}
temp <- extinction.prob(b=0.2,d=0.2, nsims=1000)
write.csv(temp,file="Tables\\RandomWalk02_quasi1000.csv")

temp <- extinction.prob(b=0.2,d=0.2, nsims=10000)
write.csv(temp,file="Tables\\RandomWalk02_quasi10000.csv")