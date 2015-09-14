#======Template: Lab 6 fitting models in R=====================================
#Using optim in R to fit a model to Antarctic blue whale abundance 
#estimates from Branch (2007). 
#Mid-year   Abundance   CV
#1981	453	0.40
#1988	559	0.47
#1998	2280	0.36
#To explore: finding MLE, exploring starting points for parameter estimates
#obtaining likelihood profiles of model parameters
#Written by Trevor A. Branch tbranch@uw.edu
#starting 6 May 2013
#====================================================================

#===PART 1=======================================
#Exponential growth function, takes parameters 
#N0 (numbers in 1975), r (rate of increase), startyear and endyear
#================================================
exp.growth <- function(N0, r, startyear=1975, endyear=2010) {

   return(N.vec)
}
exp.growth(N0=300, r=0.05, startyear=1975, endyear=2010)

#===PART 2===================================================
#Calculate the NLL of the function given the data
#Requires N0, r, startyear, endyear. The data are built-in. 
#============================================================
get.NLL <- function(N0, r, startyear=1975, endyear=2010) {
   
   return(NLL)
}
get.NLL(N0=300, r=0.05, startyear=1975, endyear=2010)


#===PART 3===================================================
#Rewrite the function so that the parameters to be estimated
#are stored in a vector.
#Requires N0, r, startyear, endyear. The data are built-in. 
#============================================================
minimize.NLL <- function(param.vector, startyear=1975, endyear=2010) {

   
   return(NLL)
}
minimize.NLL(param.vector=c(300,0.05))

#===PART 4===================================================
optim()


#============PART 5 likelihood profiling===============================================
#The main parameter of interest is r, here we will produce a likelihood profile of r
#First we need a new NLL function that has one parameter to minimize (r) and one fixed
#parameter (N0). Then we will loop over FIXED values of r and use optim() to solve for 
#the N0 value that minimizes the NLL.
#========================================================================================
Rprofile.NLL <- function(param.vector, r, startyear=1975, endyear=2010) {

   return(NLL)
}
Rprofile.NLL(param.vector=c(300), r=0.05)

optim(par=c(300),fn=Rprofile.NLL, method="Brent", r=0.05, lower=10, upper=10000)

#does the likelihood profile
Rprofile <- function(R.vec, lower=10, upper=10000) {
   

   
   
   
   plot(x=R.vec, y=saved.NLL, type="l",xaxs="i", las=1, ylab="NLL", xlab="r values", 
        lwd=2, col="blue")   
   abline(h=min(saved.NLL), lty=2, col="gray50")
   abline(h=min(saved.NLL)+1.92, lty=2, col="gray50")  #95% CI
   
   return(saved.NLL)
}
#calculate the NLL profile and plot it
R.vec <- seq(0,0.12,0.001)
x <- Rprofile(R.vec=R.vec, lower=10, upper=10000)

#calculate the likelihood profile and plot it
x.like <- vector(length=length(x))
for (i in 1:length(x)) {
   x.like[i] <- exp(min(x) - x[i])  
}
plot(x=R.vec, y=x.like, type="l",xaxs="i", yaxs="i", las=1, ylab="Scaled likelihood", xlab="r values", 
     lwd=2, col="blue")   
abline(h=exp(-1.92), lty=2, col="gray50")


#===PART 6===================================================
#testing convergence
test.optim <- function(N0.vec, r.vec, method="Nelder-Mead") {
   rownames <- paste("r=",r.vec)
   colnames <- paste("N0=",N0.vec)
   NLL.mat <- matrix(nrow=length(r.vec), ncol=length(N0.vec), 
                     dimnames=list(rownames, colnames))
   
   
   
   
   return(NLL.mat)
}
x <- test.optim(N0.vec=seq(100,1000,100), r.vec=seq(0,0.12,0.01), method="Nelder-Mead")
x <- test.optim(N0.vec=seq(100,1000,100), r.vec=seq(0,0.12,0.01), method="BFGS")
x <- test.optim(N0.vec=seq(100,1000,100), r.vec=seq(0,0.12,0.01), method="CG")
#x <- test.optim(N0.vec=seq(100,1000,100), r.vec=seq(0,0.12,0.01), method="L-BFGS-B")  #needs bounded pars
#x <- test.optim(N0.vec=seq(100,1000,100), r.vec=seq(0,0.12,0.01), method="SANN")      #takes a long time
#x <- test.optim(N0.vec=seq(100,1000,100), r.vec=seq(0,0.12,0.01), method="Brent")     #only for one param
