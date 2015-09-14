##################################################
#Individual based model in R 
#model is N(t+1) = N(t) + b*N(t) - m*N(t)
#where b is proportion of births and m is proportions of deaths
#at t = 1, N(t)=N0
#using binomial draws to speed up the individual based model
###################################################
#Written by Trevor A. Branch 2 April 2013 for FISH458
##################################################
indiv.model.binom <- function(b=0.1, m=0.1, N0=100, nyears=50) {
   N.vec <- vector(length=nyears)   #start in year 0 end in year 50
   N.vec[1] <- N0
   
   for (yr in 2:nyears) {   #loop through years
      births <- rbinom(n=1,size=N.vec[yr-1],prob=b)   #random draw from binomial (n=1 means one draw) 
      deaths <- rbinom(n=1,size=N.vec[yr-1],prob=m)
      N.vec[yr] <- N.vec[yr-1]+births-deaths
   }
   return(N.vec)
}
#plot one trajectory
Abundance <- indiv.model.binom(b=0.15,m=0.13)
plot(Abundance, ylim=c(0,1.2*max(Abundance)), type="l", yaxs="i")    #ylim, type, yaxs makes plots nicer

#plot ten trajectories
for (j in 1:10) {
   Abundance <- indiv.model.binom(b=0.15,m=0.13, N0=100, nyears=100)
   if (j>1) { 
      par(new=T) 
   }
   plot(Abundance, ylim=c(0,400), type="l", yaxs="i")    #ylim, type, yaxs makes plots nicer
}