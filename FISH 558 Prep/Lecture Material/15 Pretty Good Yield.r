#############################################################
#Model for Pretty Good Yield after the paper: 
#Hilborn, R. 2010. Pretty Good Yield and exploited fishes. 
#Mar. Pol. 34: 193-196.
#Builds an age-structured model and explores how yield changes
#relative to harvest rate for a variety of model assumptions.
#Default parameter value from Atlantic halibut from FishBase
#Many of the model options were first calculated by 
#Steve Hare, srhare@gmail.com when he taught FISH458 in 2011
#but the model has been re-coded here completely. 
#Trevor A. Branch tbranch@uw.edu
#Created 2 May 2013
#############################################################

#++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++model+++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++
calc.N.at.age <- function(R, Surv, n, u, vfull)  {
   N.vec <- vector(length=n)  #numbers at each age
   vuln <- vector(length=n)
   
   #vulnerability at age
   vuln[vfull:n] <- 1   #all other values by default are set to 0 in R
   
   N.vec[1] <- R
   for(i in 1:(n-2)) {
      N.vec[i+1] <- N.vec[i]*Surv*(1-vuln[i]*u)  
   }
   N.vec[n] <- Surv*(1-vuln[n]*u) / (1-(Surv*(1-vuln[n]*u))) * N.vec[n-1]
   
   return(N.vec)
}


calc.yield <- function(R0, Surv, vfull, am, u, alpha, beta, Linfinity, k, t0, n, pfemale, h)  {
   #calculate the number of fish in an unfished population at each age, at equilibrium
   N.unfished <- calc.N.at.age(R=R0, Surv=Surv, n=n, u=0, vfull=vfull)
   
   #length at age vector
   lengths <- vector(length=n)   #lengths in cm 
   for (i in 1:n) {
      lengths[i] <- Linfinity*(1-exp(-1*k*(i-t0)))  
   }
   
   #weight at age vector
   weights <- vector(length=n)    #weight in kg
   for (i in 1:n) {
       weights[i] <- alpha*(lengths[i]^beta)  
   }
   
   #maturity at age vector
   maturity <- vector(length=n)
   maturity[am:n] <- 1   #all other values by default are set to 0 in R

   #vulnerability at age vector
   vuln <- vector(length=n)
   vuln[vfull:n] <- 1   #all other values by default are set to 0 in R
   
   #calculate the spawning biomass in an unfished population, used as a parameter
   SSB0 <- 0  #measured in kg
   for(i in 1:n) {
      SSB0 <- SSB0 + maturity[i]*pfemale*weights[i]*N.unfished[i]
   }
   
   #calculate the parameters of the Beverton-Holt stock-recruit relation
   BH.alpha <- SSB0/R0*(1-(h-0.2)/(0.8*h))
   BH.beta <- (h-0.2)/(0.8*h*R0)
   
   #SBPR(u) spawning biomass per recruit as a function of u
   #and YPR(u) the yield per recruit as a function of u
   N.per.recruit <- calc.N.at.age(R=1, Surv=Surv, n=n, u=u, vfull=vfull)
   SBPR <- 0
   YPR <- 0
   for (i in 1:n) {
      SBPR <- SBPR + maturity[i]*pfemale*weights[i]*N.per.recruit[i]
      YPR <- YPR + vuln[i]*u*weights[i]*N.per.recruit[i]
   }
   
   #calculate equilibrium number of recruits
   eqm.recruits <- (SBPR-BH.alpha)/(BH.beta*SBPR)  #at equilibrium under this level of u
   if (eqm.recruits <0 ) { 
      eqm.recruits <- 0 
   }
   yield <- YPR*eqm.recruits
   SSB <- SBPR*eqm.recruits
   
   return(list(SSB0=SSB0, yield=yield, SSB=SSB, eqm.recruits=eqm.recruits))
}
###=====parameter values (default)
R0 <- 1000000         #recruitment in the unfished population
Surv <- 0.8           #natural survival
vfull <- 6            #age at full vulnerability to fishing (knife-edge)
am <- 6               #age at full maturity (assumed knife-edge)
u <- 0                #exploitation rate
alpha <- 0.00002760   #length-weight alpha
beta <- 2.953         #length-weight beta
Linfinity <- 220      #von Bertalanffy Linfinity
k <- 0.1              #von Bertalanffy K
t0 <- -0.1            #von Bertalanffy t0
n <- 30               #plus group age; max age is 50
pfemale <- 0.5        #proportion of population that is female
h <- 0.8              #Beverton-Holt steepness must be (0.2, 1.0)

#function calls
calc.yield(R0=R0, Surv=Surv, vfull=vfull, am=am, u=u, alpha=alpha, 
           beta=beta, Linfinity=Linfinity, k=k, t0=t0, n=n, pfemale=pfemale, h=h)

#++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++plotting of this function+++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++

###===============================
#Effect of steepness on yield
harvest.steepness <- function() {
   R0 <- 1000000         #recruitment in the unfished population
   Surv <- 0.8           #natural survival
   vfull <- 6            #age at full vulnerability to fishing (knife-edge)
   am <- 6               #age at full maturity (assumed knife-edge)
   u <- 0                #exploitation rate
   alpha <- 0.00002760   #length-weight alpha
   beta <- 2.953         #length-weight beta
   Linfinity <- 220      #von Bertalanffy Linfinity
   k <- 0.1              #von Bertalanffy K
   t0 <- -0.1            #von Bertalanffy t0
   n <- 30               #plus group age; max age is 50
   pfemale <- 0.5        #proportion of population that is female
   h <- 0.8              #Beverton-Holt steepness must be (0.2, 1.0)
   
   h.vec <- seq(0.3,0.9,0.1)
   u.vec <- seq(0,1,0.001)
   results <- matrix(nrow=length(h.vec), ncol=length(u.vec))
   for (i in 1:length(h.vec)) {
      for (j in 1:length(u.vec)) {
         results[i,j] <- calc.yield(R0=R0, Surv=Surv, vfull=vfull, am=am, u=u.vec[j], alpha=alpha, 
                 beta=beta, Linfinity=Linfinity, k=k, t0=t0, n=n, pfemale=pfemale, h=h.vec[i])$yield
      }
   }
   bluepal <- colorRampPalette(c("white","blue"))
   bluecols <- bluepal(n=length(h.vec)+2)[-c(1:2)]  #create colors and delete the last two
   plot(x=u.vec, y=results[1,]/1000000, type="l", xlab="Harvest rate (u)", ylab="Yield (1000 t)", lwd=2,las=1,
        xaxs="i", yaxs="i", ylim=c(0,1.05*max(results)/1000000), col=bluecols[1])
   for (i in 2:length(h.vec)) {
      lines(x=u.vec, y=results[i,]/1000000, lwd=2, col=bluecols[i])
   }
   
   return(results)
} 
pdf("Figs\\Harvest vs steepness.pdf",width=8, height=5)
x <- harvest.steepness()   
dev.off()


###===============================
#Effect of steepness on yield and SSB/SSB0
depletion.yield <- function() {
   R0 <- 1000000         #recruitment in the unfished population
   Surv <- 0.8           #natural survival
   vfull <- 6            #age at full vulnerability to fishing (knife-edge)
   am <- 6               #age at full maturity (assumed knife-edge)
   u <- 0                #exploitation rate
   alpha <- 0.00002760   #length-weight alpha
   beta <- 2.953         #length-weight beta
   Linfinity <- 220      #von Bertalanffy Linfinity
   k <- 0.1              #von Bertalanffy K
   t0 <- -0.1            #von Bertalanffy t0
   n <- 30               #plus group age; max age is 50
   pfemale <- 0.5        #proportion of population that is female
   h <- 0.8              #Beverton-Holt steepness must be (0.2, 1.0)
   
   h.vec <- seq(0.3,0.9,0.1)
   u.vec <- seq(0,1,0.001)
   X.mat <- matrix(nrow=length(h.vec), ncol=length(u.vec))
   Y.mat <- matrix(nrow=length(h.vec), ncol=length(u.vec))
   for (i in 1:length(h.vec)) {
      for (j in 1:length(u.vec)) {
         temp <- calc.yield(R0=R0, Surv=Surv, vfull=vfull, am=am, u=u.vec[j], alpha=alpha, 
                                    beta=beta, Linfinity=Linfinity, k=k, t0=t0, n=n, pfemale=pfemale, h=h.vec[i])
         X.mat[i,j] <- temp$SSB/temp$SSB0 
         Y.mat[i,j] <- temp$yield
      }
   }
   bluepal <- colorRampPalette(c("white","red"))
   bluecols <- bluepal(n=length(h.vec)+2)[-c(1:2)]  #create colors and delete the last two
   plot(x=X.mat[1,], y=Y.mat[1,]/1000000, type="l", xlab="SSB/SSB0", ylab="Yield (1000 t)", lwd=2,las=1,
        xaxs="i", yaxs="i", ylim=c(0,1.05*max(Y.mat)/1000000), col=bluecols[1])
   for (i in 2:length(h.vec)) {
      lines(x=X.mat[i,], y=Y.mat[i,]/1000000, lwd=2, col=bluecols[i])
   }
   
   return(results)
} 
pdf("Figs\\Depletion vs yield.pdf",width=8, height=5)
x <- depletion.yield()   
dev.off()

###===============================
#Effect of steepness on SSB/SSB0 vs harvest rate
depletion.harvestrate <- function() {
   R0 <- 1000000         #recruitment in the unfished population
   Surv <- 0.8           #natural survival
   vfull <- 6            #age at full vulnerability to fishing (knife-edge)
   am <- 6               #age at full maturity (assumed knife-edge)
   u <- 0                #exploitation rate
   alpha <- 0.00002760   #length-weight alpha
   beta <- 2.953         #length-weight beta
   Linfinity <- 220      #von Bertalanffy Linfinity
   k <- 0.1              #von Bertalanffy K
   t0 <- -0.1            #von Bertalanffy t0
   n <- 30               #plus group age; max age is 50
   pfemale <- 0.5        #proportion of population that is female
   h <- 0.8              #Beverton-Holt steepness must be (0.2, 1.0)
   
   h.vec <- seq(0.3,0.9,0.1)
   u.vec <- seq(0,1,0.001)
   Y.mat <- matrix(nrow=length(h.vec), ncol=length(u.vec))
   for (i in 1:length(h.vec)) {
      for (j in 1:length(u.vec)) {
         temp <- calc.yield(R0=R0, Surv=Surv, vfull=vfull, am=am, u=u.vec[j], alpha=alpha, 
                            beta=beta, Linfinity=Linfinity, k=k, t0=t0, n=n, pfemale=pfemale, h=h.vec[i])
         Y.mat[i,j] <- temp$SSB/temp$SSB0
      }
   }
   bluepal <- colorRampPalette(c("white","darkmagenta"))
   bluecols <- bluepal(n=length(h.vec)+2)[-c(1:2)]  #create colors and delete the last two
   print(length(u.vec))
   print(dim(Y.mat))
   plot(x=u.vec, y=Y.mat[1,], type="l", xlab="Harvest rate (u)", ylab="SSB/SSB0", lwd=2,las=1,
        xaxs="i", yaxs="i", ylim=c(0,1), col=bluecols[1])
   for (i in 2:length(h.vec)) {
      lines(x=u.vec, y=Y.mat[i,], lwd=2, col=bluecols[i])
   }
   
   return(results)
} 
pdf("Figs\\Harvest rate vs depletion.pdf",width=8, height=5)
x <- depletion.harvestrate()   
dev.off()


###===============================
#Rescale yield to max within each steepness, plot yield vs SSB/SSB0 for different steepness
rescaled.depletion.yield <- function() {
   R0 <- 1000000         #recruitment in the unfished population
   Surv <- 0.8           #natural survival
   vfull <- 6            #age at full vulnerability to fishing (knife-edge)
   am <- 6               #age at full maturity (assumed knife-edge)
   u <- 0                #exploitation rate
   alpha <- 0.00002760   #length-weight alpha
   beta <- 2.953         #length-weight beta
   Linfinity <- 220      #von Bertalanffy Linfinity
   k <- 0.1              #von Bertalanffy K
   t0 <- -0.1            #von Bertalanffy t0
   n <- 30               #plus group age; max age is 50
   pfemale <- 0.5        #proportion of population that is female
   h <- 0.8              #Beverton-Holt steepness must be (0.2, 1.0)
   
   h.vec <- seq(0.3,0.9,0.1)
   u.vec <- seq(0,1,0.001)
   X.mat <- matrix(nrow=length(h.vec), ncol=length(u.vec))
   Y.mat <- matrix(nrow=length(h.vec), ncol=length(u.vec))
   for (i in 1:length(h.vec)) {
      for (j in 1:length(u.vec)) {
         temp <- calc.yield(R0=R0, Surv=Surv, vfull=vfull, am=am, u=u.vec[j], alpha=alpha, 
                            beta=beta, Linfinity=Linfinity, k=k, t0=t0, n=n, pfemale=pfemale, h=h.vec[i])
         X.mat[i,j] <- temp$SSB/temp$SSB0 
         Y.mat[i,j] <- temp$yield
      }
   }
   bluepal <- colorRampPalette(c("white","darkgreen"))
   bluecols <- bluepal(n=length(h.vec)+2)[-c(1:2)]  #create colors and delete the last two
   plot(x=X.mat[1,], y=Y.mat[1,]/max(Y.mat[1,]), type="l", xlab="SSB/SSB0", ylab="Yield (1000 t)", lwd=2,las=1,
        xaxs="i", yaxs="i", ylim=c(0,1.05), col=bluecols[1])
   for (i in 2:length(h.vec)) {
      lines(x=X.mat[i,], y=Y.mat[i,]/max(Y.mat[i,]), lwd=2, col=bluecols[i])
   }
   abline(h=0.8,lty=2,lwd=1.7,col="black")
   
   return(results)
} 
pdf("Figs\\Rescaled depletion vs yield.pdf",width=8, height=5)
x <- rescaled.depletion.yield()   
dev.off()


###===============================
#Region of PGY plotted on steepness vs. SSB/SSB0
#NOT DONE AT ALL RIGHT
#NEED TO do differently: find for each steepness the u resulting in 80% max yield
#then calculate the SSB/SSB0 for those values of u, and plot those points 
PGY.steepness.depletion <- function() {
   R0 <- 1000000         #recruitment in the unfished population
   Surv <- 0.8           #natural survival
   vfull <- 6            #age at full vulnerability to fishing (knife-edge)
   am <- 6               #age at full maturity (assumed knife-edge)
   u <- 0                #exploitation rate
   alpha <- 0.00002760   #length-weight alpha
   beta <- 2.953         #length-weight beta
   Linfinity <- 220      #von Bertalanffy Linfinity
   k <- 0.1              #von Bertalanffy K
   t0 <- -0.1            #von Bertalanffy t0
   n <- 30               #plus group age; max age is 50
   pfemale <- 0.5        #proportion of population that is female
   h <- 0.8              #Beverton-Holt steepness must be (0.2, 1.0)
   
   h.vec <- seq(0.3,0.9,0.001)
   u.vec <- seq(0,1,0.01)
   X.mat <- matrix(nrow=length(h.vec), ncol=length(u.vec))
   Y.mat <- matrix(nrow=length(h.vec), ncol=length(u.vec))
   Y.scaled <- Y.mat
   plot.size <- Y.scaled
   for (i in 1:length(h.vec)) {
      for (j in 1:length(u.vec)) {
         temp <- calc.yield(R0=R0, Surv=Surv, vfull=vfull, am=am, u=u.vec[j], alpha=alpha, 
                            beta=beta, Linfinity=Linfinity, k=k, t0=t0, n=n, pfemale=pfemale, h=h.vec[i])
         X.mat[i,j] <- temp$SSB/temp$SSB0 
         Y.mat[i,j] <- temp$yield
      }
      Y.scaled[i,] <- Y.mat[i,]/max(Y.mat[i,])   #relative to maximum yield at that steepness
   }
   plot.size[] <- 0
   plot.size[Y.scaled>=0.8] <- 1
   bluepal <- colorRampPalette(c("white","darkgreen"))
   bluecols <- bluepal(n=length(h.vec)+2)[-c(1:2)]  #create colors and delete the last two
   plot(x=rep(h.vec[1],length(u.vec)), y=X.mat[1,], type="p", xlab="Steepness (h)", ylab="Depletion (SSB/SSB0)", lwd=2,las=1,
        xaxs="i", yaxs="i", ylim=c(0,1.05), col=bluecols[1], xlim=c(0.3,0.9))
   for (i in 2:length(h.vec)) {
      points(x=rep(h.vec[i],length(u.vec)), y=X.mat[i,], col=bluecols[i])
   }
   
   return(Y.scaled)
} 
pdf("Figs\\PGY steepness depletion.pdf",width=8, height=5)
x <- PGY.steepness.depletion()   
dev.off()
   
   