
# The Goal ----------------------------------------------------------------
# The purpose of this model is to use MCMC and the Leplace's Demon function to
# fit a hierarchecal model of weight-at-length to 20 populations of some random critter
library(LaplacesDemon)
library(ggplot2)
library(dplyr)

rm(list = ls())

TheData <- as.data.frame(matrix(scan("Data.dat"),ncol=3,byrow=3))

colnames(TheData) <- c('pop','length','weight')

ggplot(TheData,aes( x = length, y = weight)) + geom_point()

testfit <- function(par,Data)
{
  predweight <- par[1] * Data$length ^ par[2]

  nll <- -dnorm(log(Data$weight), log(predweight),par[3], log = T)

  return(nll)
}

mon.names <- c("LP")
parm.names <- as.parm.names(list(a=rep(1,20),b=rep(1,20),mu.a = 1, mu.b = 1, sigma.a = 1, sigma.b = 1, sigma = 1 ))

# parm <- c(rep(guess$par[1],20), rep(guess$par[2],20), rep(1,5))
parm <- c(rep(log(0.01), 20), rep(log(3),20), log(0.01),log(3), log(0.05), log(0.1), log(1))

# parm.names <- as.parm.names(list(mean.alpha=log(0.01),sigma.alpha=log(0.05),mean.beta=log(3),sigma.beta=log(0.05),alphas=rep(log(0.01),Npop),betas=rep(log(3),Npop),log.sigma=0))
parm.names <- as.parm.names(list(a=rep(log(0.01), 20),b=rep(log(3),20),mu.a = log(0.01), mu.b = log(3), sigma.a = log(0.05), sigma.b = log(0.1), sigma = 0 ))



MyData <- list(N=dim(TheData)[1] * 10,Data = TheData,
               mon.names=mon.names, parm.names=parm.names) #Mandatory object for L's D

# N <- dim(TheData)[1]                              # Number of data points
# J <- 45
#

Model <- function(parm,Data)
{
  #Put in parameters
  dat <- Data$Data

#   mean.alpha <- parm[1]
#   sigma.alpha <- exp(parm[2])
#   mean.beta <- parm[3]
#   sigma.beta <- exp(parm[4])
#   alphas <- parm[5:(4+Npop)]
#   betas <- parm[(5+Npop):(4+2*Npop)]
#   sigma <- exp(parm[5+2*Npop])
#   Alphas <- exp(alphas)
#   Betas  <- exp(betas)
#
  a <- exp(parm[1:20])

  b <- exp(parm[21:40])

  mu.a <- (parm[41])

  mu.b <- (parm[42])

  sigma.a <- exp(parm[43])

  sigma.b <- exp(parm[44])

  sigma <- exp(parm[45])

  a.frame <- data_frame(pop = 1:length(a),a = a)

  b.frame <- data_frame(pop = 1:length(b),b = b)
  dat <- dat %>%
    full_join(a.frame, by = 'pop') %>%
    full_join(b.frame, by = 'pop')

  # Specify the prior

  mu.a.prior <- dnormv(mu.a,0,1000,log=T)               # Uninformative

  mu.b.prior <- dnormv(mu.b,0,1000,log=T)                 # Uninformative

  sigma.a.prior <- dhalfcauchy(sigma.a, 25, log=T)   # One of many priors

  sigma.b.prior <- dhalfcauchy(sigma.b, 25, log=T)   # One of many priors

  sigma.prior <- dhalfcauchy(sigma, 25, log=T)   # One of many priors

  a.prior <- dnormv(log(a),mu.a,sigma.a,log=T)           # One of many priors

  b.prior <- dnormv(log(b),mu.b,sigma.b, log = T)     # One of many priors

#   alphas.prior <- dnormv(alphas,mean.alpha,sigma.alpha,log=T)       # Informative
#   betas.prior <- dnormv(betas,mean.beta,sigma.beta,log=T)

  # Compute the log-likelihood

  pred.weight <- with(dat,a * (length)^b) #Predicted
#   pred.weight <- with(dat,log(a) + b*log(length)) #Predicted

  LL <- sum(dnorm(log(dat$weight), log(pred.weight), sigma, log=TRUE))     # Log-likelihood

  # Log-posterior
  LP <- LL + sum(a.prior) + sum(b.prior) +  mu.a.prior + mu.b.prior +  sigma.a.prior + sigma.b.prior+ sigma.prior
  # Return key stuff
  Modelout <- list(LP=LP,Dev=-2*LL, Monitor= c(LP), yhat=(rnorm(length(pred.weight), pred.weight, sigma)),parm=parm)
  return(Modelout)
}

mine <- Model(parm = parm, Data = MyData)

Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values = parm,Covar=NULL,
                     mon.names = mon.names,Iterations=2500000, Status=10000,
                     Thinning=20000,Algorithm="HARM", Specs=NULL)

# Initial.Values <- jitter(as.initial.values(Fit))
#
# huh <- Model(parm = Initial.Values, Data = MyData)
#
# huh$yhat
#
# plot(TheData$weight,huh$yhat)
# abline(a=0, b=1)
#
# Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#                      Covar=NULL, Iterations=3e6, Status=65250, Thinning=1000,
#                      Algorithm="HARM")
#
# huh <- Model(parm = Initial.Values, Data = MyData)
#
# huh$yhat
#
# plot(TheData$weight,huh$yhat)
# abline(a=0, b=1)
#

Consort(Fit)
print(Fit)
plot(Fit,BurnIn=20000,MyData,PDF=F,Parms=NULL)
#
# # Plots and diagnostics
# plot(Fit,BurnIn=100,MyData,PDF=FALSE,Parms=NULL)
# plot(Fit.hpc,BurnIn=100,MyData,PDF=FALSE,Parms=NULL)
# caterpillar.plot(Fit.hpc,Parms=1:3)
# Pred <- predict(Fit,Model,MyData)
# summary(Pred,Discrep="Chi-Square")
# plot(Pred,Style="Density",Rows=1:9)
# plot(Pred,Style="Fitted")
# plot(Pred,Style="Jarque-Bera")
# plot(Pred,Style="Residuals")



