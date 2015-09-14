##################################################
#Individual based model in R 
#model is N(t+1) = N(t) + b*N(t) - m*N(t)
#where b is proportion of births and m is proportions of deaths
#at t = 1, N(t)=N0
#this is the long-winded way of going through every single individual
#and randomly deciding if it gives birth and/or dies
##################################################
#Written by Trevor A. Branch 2 April 2013 for FISH458
##################################################
indiv.model <- function(b=0.1, m=0.1, N0=100, nyears=50) {
   N.vec <- vector(length=nyears)   #start in year 0 end in year 50
   N.vec[1] <- N0
   
   for (yr in 2:nyears) {   #loop through years
      births <- 0
      deaths <- 0
      for (i in 1:N.vec[yr-1])  {
         if(runif(n=1) < b) {   #if random number smaller than birth rate then new birth
            births <- births+1
         }   
         if(runif(n=1) < m) {   #if random number smaller than birth rate then new birth
            deaths <- deaths+1
         }   
      }
      N.vec[yr] <- N.vec[yr-1]+births-deaths
   }
   return(N.vec)
}
#plot one trajectory
Abundance <- indiv.model(b=0.15,m=0.13)
plot(Abundance, ylim=c(0,1.2*max(Abundance)), type="l", yaxs="i")    #ylim, type, yaxs makes plots nicer

#plot ten trajectories
for (j in 1:10) {
   Abundance <- indiv.model(b=0.15,m=0.13, N0=100, nyears=100)
   if (j>1) { 
      par(new=T) 
   }
   plot(Abundance, ylim=c(0,400), type="l", yaxs="i")    #ylim, type, yaxs makes plots nicer
}