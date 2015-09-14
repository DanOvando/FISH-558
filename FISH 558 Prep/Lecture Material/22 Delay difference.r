###########################################################################
#Exploration of different types of delay-difference models
#using the manipulate() function to explore the options.
#Written by Trevor A. Branch tbranch@uw.edu
#Starting 21 May 2013
###########################################################################

######Delay-logistic  does not seem to work as promised####################
#r=rate of increase, K=carrying capacity, S=survival, 
#L=lag from birth to recruitment
#scenario: equilibrium at N0 (perhaps from harvest) from year 1 to year L
#then no harvest and allowed to return to K. 
#Note 1: S=1 unless there is additional mortality (e.g. fishing)
#Note 2: when r*L > 1 there will be a limit cycle, even with smaller 
#        r there will be damped oscillations
delay.logistic <- function(N0, r, K, S, L, nyears) {
   N <- vector(length=nyears+L+1)
   N[1:(L+1)] <- N0   #first L years it is stable at N0
   
   for (i in (L+1):(nyears+L)) {
      N[i+1] <- S*N[i] + r*N[i-L]*(1-N[i-L]/K)
   }
   plot(x=1:(nyears+L+1), y=N, type="l",lwd=2,col="blue", ylim=c(0,1.2*K))   
   
   #compare with regular logistic
   N <- vector(length=nyears+L+1)
   N[1:(L+1)] <- N0   #first L years it is stable at N0
   
   for (i in (L+1):(nyears+L)) {
      N[i+1] <- N[i] + r*N[i]*(1-N[i]/K)
   }
   lines(x=1:(nyears+L+1), y=N, lwd=2,col="green")   
}
par(mfrow=c(3,1), mar=c(2,2,2,2), oma=c(1,1,0,0))
delay.logistic(N0=100, r=0.05, K=1000, S=1, L=5, nyears=150)
delay.logistic(N0=100, r=0.2, K=1000, S=1, L=5, nyears=150)
delay.logistic(N0=100, r=0.3, K=1000, S=1, L=5, nyears=150)



###############Deriso-Schnute delay-difference model##########################
#s = total survival, S = natural survival, C = catches, rho = growth coefficient,
#WL = weight at age L, WLm1 = weight at age L-1, 
#Rt = recruitment in year t, Rtplus1 = recruitment in time t+1
#R0 = unfished equilibrium recruitment
#B0 = function of R0 and the other parameters
##############################################################################
Deriso.Schnute <- function() {
   #For those interested: how would you implement the Deriso-Schnute delay-difference? 
   
}