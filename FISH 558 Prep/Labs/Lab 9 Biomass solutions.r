#======FISH458 Lab 9===================================================================
#Understanding simulation-estimation methods. Set up an operating model that generates
#simulated data and true biomass estimates. Then apply a management rule to the 
#simulated data and test whether the management rule properly stabilizes the abundance 
#and produces reasonable catches. 
#Written by Trevor A. Branch, tbranch@uw.edu
#starting 27 May 2013.
#With help from earlier code by Cole Monnahan.
#======================================================================================


#=======================================================================================
#Runs the model to year 21 and returns the sum of the squares of the difference
#between B21/K and the required tuning level (which is 0.2*K). In other words, 
#if B21 = 0.2K then it will return zero. Used in optim to find harvest rate resulting
#in the required tuning level. 
#=======================================================================================
find.u <- function(u, r=0.3, K=1000, sigma.w=0.1, tuning.level=0.2, data.years=20, 
                   random.seed) {
   #generate 20 process errors (will be fixed for a given random.seed)
   set.seed(random.seed)
   wt <- rnorm(n=data.years,mean=0,sd=sigma.w)
   Biomass <- K
   for (i in 1:data.years) {
      Biomass= max(0,(Biomass+r*Biomass*(1-Biomass/K) - u*Biomass)*exp(wt[i]-0.5*sigma.w^2))
   }
   return((Biomass/K-tuning.level)^2)
}

#====================================================================================
#Starting biomass in the first 20 years, including finding the correct
operating.model <- function(r=0.3, K=1000, q=0.001, sigma.w=0.1, sigma.obs=0.3, tuning.level=0.2, 
                            data.years=20, nyears=35, random.seed) {
   #Find the harvest rate resulting in B/K = tuning.level in data.years+1
   #Note that by adding $par[1] we only store the estimate for u and don't do any convergence checking
   tuning.u <- optim(par=0.1,method="Brent", fn=find.u, lower=0.01, upper=0.99, 
         r=r, K=K, sigma.w=sigma.w, tuning.level=tuning.level,
         data.years=data.years,random.seed=random.seed)$par[1]

   #-----Simulate true state of nature
   #Vectors to store biomass, catch, index values in each year
   Bvec <- vector(length=nyears)
   Cvec <- vector(length=nyears)
   
   #generate 20 process errors (will be fixed for a given random.seed)
   set.seed(random.seed)
   wt <- rnorm(n=nyears,mean=0,sd=sigma.w)

   Bvec[1] <- K
   for (i in 1:data.years) {
      Cvec[i] <- tuning.u*Bvec[i]
      Bvec[i+1]= max(0,(Bvec[i]+r*Bvec[i]*(1-Bvec[i]/K) - Cvec[i])*exp(wt[i]-0.5*sigma.w^2))
   }
   
   #-----Generate simulated data from truth plus noise
   obsB <- vector(length=nyears)  #Biomass estimates ("observed biomass")
   obsI <- vector(length=nyears)  #Index estimates ("observed index values")
   obs.B.err <- rnorm(n=nyears,mean=0,sd=sigma.obs)  #observation error on B 1-21
   obs.I.err <- rnorm(n=nyears,mean=0,sd=sigma.obs)  #observation error on q 1-21
   
   for (i in 1:(data.years+1)) { 
      obsB[i] <- Bvec[i]*exp(obs.B.err[i]-0.5*sigma.obs^2)     
      obsI[i] <- Bvec[i]*q*exp(obs.I.err[i]-0.5*sigma.obs^2)
   }
   
   #----return the true values B and C, and the simulated data obsB and obsI
   return(list(B=Bvec, C=Cvec, obsB=obsB, obsI=obsI, wt=wt, 
               obs.B.err=obs.B.err, obs.I.err=obs.I.err))   #notice I am naming the elements of the list
}
temp <- operating.model(random.seed=5)
temp$obsB
temp[[3]]   #this is the same as saying $obsB
temp$obsI
temp$C
temp$wt     

#====PART 1=====================================================================================
#This implements the floor-rate harvest control policy, with 
#alpha=minimum biomass for an open fishery, and beta=harvest rate for 
#biomass above alpha
#===============================================================================================
project.biomass <- function(r=0.3, K=1000, q=0.001, sigma.w=0.1, sigma.obs=0.3, tuning.level=0.2, 
                data.years=20, nyears=35, random.seed, alpha=0, beta=0.15) {
   #first calculate the data for the first 20/21 years to get the model tuned 
   #note that the B, I, Bobs etc are already returned as vectors of length nyears with space for projections
   temp <- operating.model(r=r, K=K, q=q, sigma.w=sigma.w, sigma.obs=sigma.obs, tuning.level=tuning.level, 
                           data.years=data.years, nyears=nyears, random.seed=random.seed)
   Bvec <- temp$B                #true biomass values
   Cvec <- temp$C                #catches in each year
   obsB <- temp$obsB             #observed biomass estimates
   wt <- temp$wt                 #vector of process errors
   obs.B.err <- temp$obs.B.err   #vector of observation errors for biomass estimates
   
   #now loop over future years, calculate catch, Bvec, obsB
   for (i in (data.years+1):(nyears-1))  {
      Cvec[i] <- max(0,beta*(obsB[i]-alpha))
      Bvec[i+1]= max(0,(Bvec[i]+r*Bvec[i]*(1-Bvec[i]/K) - Cvec[i])*exp(wt[i]-0.5*sigma.w^2))
      obsB[i+1] <- Bvec[i+1]*exp(obs.B.err[i+1]-0.5*sigma.obs^2)
   }
   #calculate catch in the final year
   Cvec[nyears] <- max(0,beta*(obsB[nyears]-alpha))
   
   #calculate average catch during management years
   Cmean <- mean(Cvec[(data.years+1):nyears])
   
   #calculate CV of catch during management years
   C.CV <- sd(Cvec[(data.years+1):nyears])/mean(Cvec[(data.years+1):nyears])

   #calculate mean of *true* biomass during management years
   Bmean <- mean(Bvec[(data.years+1):nyears])
   
   #calculate minimum of true biomass during management years
   Bmin <- min(Bvec[(data.years+1):nyears])
   
   #return a list of the results
   return(list(B=Bvec, C=Cvec, obsB=obsB, Cmean=Cmean, C.CV=C.CV, 
               Bmean=Bmean, Bmin=Bmin))
}
x <- project.biomass(random.seed=1, alpha=0, beta=0.15)
x

#====PART 2=====================================================================================
#Calls the function n times to repeat the calculations and plot histograms of the key 
#results: mean C, mean B, CV of catch.
#===============================================================================================
evaluate.rule <- function(alpha=0, beta=0.15, nruns=100, r=0.3, K=1000, q=0.001, sigma.w=0.1, sigma.obs=0.3, tuning.level=0.2, 
                            data.years=20, nyears=35) {
   #create vectors to store Cmean, C.CV, Bmean, Bmin
   Cmean <- vector(length=nruns)
   C.CV <- vector(length=nruns)
   Bmean <- vector(length=nruns)
   Bmin <- vector(length=nruns)
   
   #loop through nruns storing the results at each step
   for (i in 1:nruns) {
      temp <- project.biomass(r=r, K=K, q=q, sigma.w=sigma.w, sigma.obs=sigma.obs, tuning.level=tuning.level, 
                      data.years=data.years, nyears=nyears, random.seed=i, alpha=alpha, beta=beta)
      Cmean[i] <- temp$Cmean
      C.CV[i] <- temp$C.CV
      Bmean[i] <- temp$Bmean
      Bmin[i] <- temp$Bmin
   }
   
   #plot the distribution of each key parameter using hist
   par(mfcol=c(2,2), oma=c(0,0,0,0),mar=c(4,4,1,1), yaxs="i", xaxs="i")
   hist(Cmean, breaks=seq(0,1.1*max(Cmean),2), main="", col="gray")
   hist(C.CV, breaks=seq(0,1.1*max(C.CV),0.02), main="", col="gray")
   hist(Bmean, breaks=seq(0,K,50), main="", col="gray")
   hist(Bmin, breaks=seq(0,K,50), main="", col="gray")
}
evaluate.rule(nruns=100, alpha=0,beta=0.15)

#explore the results using manipulate()
require(manipulate)
manipulate(evaluate.rule(alpha, beta, nruns, r=0.3, K=1000, q=0.001, sigma.w=0.1, sigma.obs=0.3, tuning.level=0.2, 
                         data.years=20, nyears=35), 
           alpha=slider(0,400, initial=0),
           beta=slider(0.01,0.6, initial=0.15),
           nruns=slider(100,2000, initial=100))

