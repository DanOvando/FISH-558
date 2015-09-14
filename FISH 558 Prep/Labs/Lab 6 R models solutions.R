#======Solutions: Lab 6 fitting models in R=====================================
#Using optim in R to fit a model to Antarctic blue whale abundance 
#estimates from Branch (2007). 
#Mid-year   Abundance	CV
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
   nyears <- endyear-startyear+1
   N.vec <- vector(length=nyears)  
   N.vec[1] <- N0
   for (i in 1:(nyears-1)) {
      N.vec[i+1] <- (1+r)*N.vec[i]  
   }
   return(N.vec)
}
exp.growth(N0=300, r=0.05, startyear=1975, endyear=2010)

#===PART 2===================================================
#Calculate the NLL of the function given the data
#Requires N0, r, startyear, endyear. The data are built-in. 
#============================================================
get.NLL <- function(N0, r, startyear=1975, endyear=2010) {
   nyears <- endyear-startyear+1
   N.vec <- vector(length=nyears)  
   N.vec[1] <- N0
   for (i in 1:(nyears-1)) {
      N.vec[i+1] <- (1+r)*N.vec[i]  
   }
   
   N1981 <- N.vec[1981-startyear+1] 
   N1988 <- N.vec[1988-startyear+1] 
   N1998 <- N.vec[1998-startyear+1] 
   print(N1981)
   print(N1988)
   print(N1998)
   NLL <- 0
   NLL <- (log(453)- log(N1981))^2/(2*0.40^2) +
          (log(559)- log(N1988))^2/(2*0.47^2) +   
          (log(2280)-log(N1998))^2/(2*0.36^2)
   
   return(NLL)
}
get.NLL(N0=300, r=0.05, startyear=1975, endyear=2010)


#===PART 3===================================================
#Rewrite the function so that the parameters to be estimated
#are stored in a vector.
#Requires N0, r, startyear, endyear. The data are built-in. 
#============================================================
minimize.NLL <- function(param.vector, startyear=1975, endyear=2010) {
   N0 <- param.vector[1]
   r <- param.vector[2]
   nyears <- endyear-startyear+1
   N.vec <- vector(length=nyears)  
   N.vec[1] <- N0
   for (i in 1:(nyears-1)) {
      N.vec[i+1] <- (1+r)*N.vec[i]
   }
   
   N1981 <- N.vec[1981-startyear+1] 
   N1988 <- N.vec[1988-startyear+1] 
   N1998 <- N.vec[1998-startyear+1] 
   #print(N1981)
   #print(N1988)
   #print(N1998)
   NLL <- 0
   NLL <- (log(453)- log(N1981))^2/(2*0.40^2) +
      (log(559)- log(N1988))^2/(2*0.47^2) +   
      (log(2280)-log(N1998))^2/(2*0.36^2)
   
   return(NLL)
}
minimize.NLL(param.vector=c(300,0.05))

#===PART 4===================================================
optim(par=c(300,0.05),fn=minimize.NLL, method="Nelder-Mead")

optim(par=c(300,0.05),fn=minimize.NLL, method="BFGS")
optim(par=c(300,0.05),fn=minimize.NLL, method="BFGS", 
      control=list(maxit=5000))   #try increasing the number of iterations

#============PART 5 likelihood profiling===============================================
#The main parameter of interest is r, here we will produce a likelihood profile of r
#First we need a new NLL function that has one parameter to minimize (r) and one fixed
#parameter (N0). Then we will loop over FIXED values of r and use optim() to solve for 
#the N0 value that minimizes the NLL.
#========================================================================================
Rprofile.NLL <- function(param.vector, r, startyear=1975, endyear=2010) {
   N0 <- param.vector[1]
   nyears <- endyear-startyear+1
   N.vec <- vector(length=nyears)  
   N.vec[1] <- N0
   for (i in 1:(nyears-1)) {
      N.vec[i+1] <- (1+r)*N.vec[i]
   }
   
   N1981 <- N.vec[1981-startyear+1] 
   N1988 <- N.vec[1988-startyear+1] 
   N1998 <- N.vec[1998-startyear+1] 
   #print(N1981)
   #print(N1988)
   #print(N1998)
   NLL <- 0
   NLL <- (log(453)- log(N1981))^2/(2*0.40^2) +
      (log(559)- log(N1988))^2/(2*0.47^2) +   
      (log(2280)-log(N1998))^2/(2*0.36^2)
   
   return(NLL)
}
Rprofile.NLL(param.vector=c(300), r=0.05)

optim(par=c(300),fn=Rprofile.NLL, method="Brent", r=0.05, lower=10, upper=10000)

#does the likelihood profile
Rprofile <- function(R.vec, lower=10, upper=10000) {
   nR <- length(R.vec)
   saved.NLL <- vector(length=nR)
   
   for (i in 1:nR) {
      x <- optim(par=c(300),fn=Rprofile.NLL, method="Brent", 
                        r=R.vec[i], lower=10, upper=10000)
      saved.NLL[i] <- x$value
   }
   plot(x=R.vec, y=saved.NLL, type="l",xaxs="i", las=1, ylab="NLL", xlab="r values", 
        lwd=2, col="blue")   
   abline(h=min(saved.NLL), lty=2, col="gray50")
   abline(h=min(saved.NLL)+1.92, lty=2, col="gray50")  #95% CI
   
   return(saved.NLL)
}
#calculate the NLL profile and plot it
R.vec <- seq(0,0.12,0.0001)
x <- Rprofile(R.vec=R.vec, lower=10, upper=10000)

#calculate the likelihood profile and plot it
x.like <- vector(length=length(x))
for (i in 1:length(x)) {
   x.like[i] <- exp(min(x) - x[i])  
}
plot(x=R.vec, y=x.like, type="l",xaxs="i", yaxs="i", 
     las=1, ylab="Scaled likelihood", xlab="r values", 
     lwd=2, col="blue")   
abline(h=exp(-1.92), lty=2, col="gray50")


#===PART 6===================================================
#testing convergence
test.optim <- function(N0.vec, r.vec, method="Nelder-Mead") {
   rownames <- paste("r=",r.vec)
   colnames <- paste("N0=",N0.vec)
   NLL.mat <- matrix(nrow=length(r.vec), ncol=length(N0.vec), 
                     dimnames=list(rownames, colnames))
   for (i in 1:length(r.vec)) {
      for (j in 1:length(N0.vec)) {
         x <- optim(par=c(N0.vec[j], r.vec[i]), fn=minimize.NLL, 
                    method=method)
         #this option increases the number of function calls which can improve convergence 
         #performance for the BFGS method. This might take a very long time (>20 minutes)
         #x <- optim(par=c(N0.vec[j], r.vec[i]), fn=minimize.NLL, 
         #           method=method, control=list(maxit=5000))
         NLL.mat[i,j] <- x$value 
         #NLL.mat[i,j] <- x$par[1] 
      }
   }
   return(NLL.mat)
}
x <- test.optim(N0.vec=seq(100,1000,100), r.vec=seq(0,0.12,0.01), method="Nelder-Mead")
x <- test.optim(N0.vec=seq(100,1000,100), r.vec=seq(0,0.12,0.01), method="BFGS")
x <- test.optim(N0.vec=seq(100,1000,100), r.vec=seq(0,0.12,0.01), method="CG")
#x <- test.optim(N0.vec=seq(100,1000,100), r.vec=seq(0,0.12,0.01), method="L-BFGS-B")  #needs bounded pars
#x <- test.optim(N0.vec=seq(100,1000,100), r.vec=seq(0,0.12,0.01), method="SANN")      #takes a long time
#x <- test.optim(N0.vec=seq(100,1000,100), r.vec=seq(0,0.12,0.01), method="Brent")     #only for one param

