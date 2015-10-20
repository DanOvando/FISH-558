library(LaplacesDemon)
set.seed(-1801)

Nsamp <- 100
X <- rnorm(Nsamp,0,5)
X <- rbind(1,X)
Mean <- c(1,4)%*%X
Y <- rnorm(Nsamp,Mean,3)

# N <- length(Y)                              # Number of data points
# J <- 2                                      # Number of parameters
mon.names <- c("LP","sigma")                 #Parameters to monitor
parm.names <- as.parm.names(list(alpha=0,beta=0,log.sigma=0))
MyData <- list(mon.names=mon.names,parm.names=parm.names,X=X,Y=Y)

Model <- function(parm,Data)
 {
  beta <- parm[1:2]
  sigma <- exp(parm[3])

  # Specify the prior
  beta.prior <- dnormv(beta,0,1000,log=T)           # Uninformative
  sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)   # One of many priors
  mu <- beta%*%Data$X

  # Compute the log-likelihood
  LL <- sum(dnorm(Data$Y, mu, sigma, log=TRUE))     # Log-likelihood

  # Log-posterior
  LP <- LL + sum(beta.prior) + sigma.prior

  # Return key stuff
  Modelout <- list(LP=LP,Dev=-2*LL, Monitor= c(LP,sigma), yhat=rnorm(length(mu), mu, sigma),parm=parm)
  return(Modelout)
}

# Fit the model
Initial.Values <- c(rep(0,2), log(5))
Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,Covar=NULL, Iterations=150000, Status=50000, Thinning=150,Algorithm="HARM", Specs=NULL)

print(Fit)
Fit.hpc <- LaplacesDemon.hpc(Model, Data=MyData, Initial.Values,Covar=NULL, Iterations=150000, Status=50000, Thinning=150,Algorithm="HARM", Specs=NULL,Chains=4,CPUs=4,Packages=NULL,Dyn.libs=NULL)
print(Fit.hpc)
Fit2 <- LaplacesDemon(Model, Data=MyData, Initial.Values,Covar=NULL, Iterations=150000, Status=50000, Thinning=150,Algorithm="AM", Specs=list(Adaptive=50000,Periodicity=10000))
print(Fit2)


# Plots and diagnostics
plot(Fit,BurnIn=100,MyData,PDF=FALSE,Parms=NULL)
plot(Fit.hpc,BurnIn=100,MyData,PDF=FALSE,Parms=NULL)
caterpillar.plot(Fit.hpc,Parms=1:3)
AAA
Pred <- predict(Fit,Model,MyData)
summary(Pred,Discrep="Chi-Square")
plot(Pred,Style="Density",Rows=1:9)
plot(Pred,Style="Fitted")
plot(Pred,Style="Jarque-Bera")
plot(Pred,Style="Residuals")






