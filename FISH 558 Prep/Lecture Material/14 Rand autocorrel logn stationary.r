#############################################################
#Stationary catches over time with mean mu but also with 
#autocorrelation with correlation rho. Error is lognormal. 
#To test correctness, take mean(N vector)
#note: simulates nyears+10 years and then deletes the first 10 years
#plot=T draws plot, plot=F does not. 
#slightly modified from code by Michael Wilberg
#Trevor A. Branch tbranch@uw.edu
#Created 2009
#Last modified 30 April 2013:
#     changed variable names for FISH 458)
#     modified plotting code
#############################################################
st.autocorrel.catch <- function(nyears=100, mu=1000,rho=0.5,CV=0.2,seed=7,plot=T) {
   set.seed(seed)   #picks a particular sequence of random numbers
   sim.years <- nyears+10
   S <- sqrt(log(CV^2+1))  #SD parameter of lognormal distribution with desired variability
   Se <- sqrt(S^2/(1-rho^2))  #SD for first value
   eps <- rnorm(mean=0,sd=1,n=sim.years)
   X <- vector(length=sim.years)
   X[1] <- Se*eps[1]
   for (i in 1:(sim.years-1)) {
      X[i+1] <- rho*X[i]+eps[i]*S
   }
   N <- exp(X+log(mu)-0.5*Se^2)
   if (plot==T) {
      plot(N,type="l",ylim=c(0,1.05*max(N)),las=1,yaxs="i", xaxs="i",
           xlab="Year",ylab="Abundance")
      abline(h=mu,lty=2,lwd=2,col="gray50")
      abline(h=mean(N),lty=2,lwd=2,col="red")
   }
   invisible(N[-c(1:10)])
}
x <- st.autocorrel.catch(nyears=10, mu=1000, rho=0.5, CV=0.1, seed=4,plot=T)
mean(x)   #should be equal to mu for large nyears
sd(x)/mean(x)   #should be equal to CV for large nyears and small rho
cor(x[-1],x[-length(x)])   #should be equal to rho for large nyears and small CV 

library(manipulate)
manipulate(st.autocorrel.catch(nyears, mu, rho, CV, seed), 
           nyears=slider(5,1000, initial=50),
           mu=slider(100,5000, initial=1000),
           rho=slider(-0.99,0.99, initial=0.1),
           CV=slider(0,1, initial=0.1),
           seed=slider(1,50,initial=7))

#create plot comparing results for three values of CV with same seed   
some.seed <- 6
nyears <- 20
x1 <- st.autocorrel.catch(nyears=nyears, mu=1000, rho=0.5, CV=0.1, seed=some.seed, plot=F)
x2 <- st.autocorrel.catch(nyears=nyears, mu=1000, rho=0.5, CV=0.25, seed=some.seed, plot=F)
x3 <- st.autocorrel.catch(nyears=nyears, mu=1000, rho=0.5, CV=0.4, seed=some.seed, plot=F)
plot(x=1:nyears,y=x1,ylim=c(0,max(c(x1,x2,x3))), yaxs="i",xaxs="i", type="l",lwd=2, col="blue",
     xlab="Years",ylab="Abundance")        
lines(x=1:nyears,y=x2, lwd=2, col="red")
lines(x=1:nyears,y=x3, lwd=2, col="green")
     
#create plot comparing results for three values of CV with same seed   
some.seed <- 6
nyears <- 20
x1 <- st.autocorrel.catch(nyears=nyears, mu=1000, rho=0, CV=0.3, seed=some.seed, plot=F)
x2 <- st.autocorrel.catch(nyears=nyears, mu=1000, rho=0.5, CV=0.3, seed=some.seed, plot=F)
x3 <- st.autocorrel.catch(nyears=nyears, mu=1000, rho=0.9, CV=0.3, seed=some.seed, plot=F)
plot(x=1:nyears,y=x1,ylim=c(0,max(c(x1,x2,x3))), yaxs="i",xaxs="i", type="l",lwd=2, col="blue",
    xlab="Years",ylab="Abundance")        
lines(x=1:nyears,y=x2, lwd=2, col="red")
lines(x=1:nyears,y=x3, lwd=2, col="green")
           