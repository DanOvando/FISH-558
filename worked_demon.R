rm(list = ls())
library(LaplacesDemon)
GenerateData <- function()
{
  set.seed(777)
  Npop <- 20
  Nsamp <- 15
  LogA <- log(0.01) + rnorm(Npop,0,0.05)
  LogB <- log(3) + rnorm(Npop,0,0.1)
  print(LogA)
  print(LogB)
  sigma <- 0.2

  YY <- NULL
  for (Ipop in 1:Npop)
  {
    Lengths <- runif(Nsamp,0,100)
    XX <- cbind(rep(Ipop,Nsamp),Lengths,exp(LogA[Ipop])*Lengths^exp(LogB[Ipop])*exp(rnorm(Nsamp,0,sigma)))
    YY <- rbind(YY,XX)
  }
  ZZ <- cbind(LogA,LogB,exp(LogA),exp(LogB))
  write(t(YY),"G:\\Data.DAT",ncolumns=3)
  write(t(ZZ),"G:\\True.DAT",ncolumns=4)

}

Analyze <- function()
{
  Npop <- 20
  Nsamp <- 15
  TheData <- read.table("G:\\Data.DAT",header=F)
  Index <- TheData[,1]
  Lengths <- TheData[,2]
  Weights <- TheData[,3]

  mon.names <- c("LP","sigma","sigma.alpha","sigma.beta")                 #Parameters to monitor
  parm.names <- as.parm.names(list(mean.alpha=log(0.01),sigma.alpha=log(0.05),mean.beta=log(3),sigma.beta=log(0.05),alphas=rep(log(0.01),Npop),betas=rep(log(3),Npop),log.sigma=0))
  MyData <- list(N=Nsamp*Npop,mon.names=mon.names,parm.names=parm.names,Index=Index,Lengths=Lengths,Weights=Weights)

  Model <- function(parm,Data)
  {
    a <- proc.time()
    Npop <- 20
    Nsamp <- 15
    # Extract parmeters
    mean.alpha <- parm[1]
    sigma.alpha <- exp(parm[2])
    mean.beta <- parm[3]
    sigma.beta <- exp(parm[4])
    alphas <- parm[5:(4+Npop)]
    betas <- parm[(5+Npop):(4+2*Npop)]
    sigma <- exp(parm[5+2*Npop])
    Alphas <- exp(alphas)
    Betas  <- exp(betas)
    #print(alphas)
    #print(betas)
    #print(sigma)
    #print(parm)

    # log-prior
    mean.alpha.prior <- dnormv(mean.alpha,0,1000,log=T)               # Uninformative
    mean.beta.prior <- dnormv(mean.beta,0,1000,log=T)                 # Uninformative
    alphas.prior <- dnormv(alphas,mean.alpha,sigma.alpha,log=T)       # Informative
    betas.prior <- dnormv(betas,mean.beta,sigma.beta,log=T)           # Informative
    sigma.prior1 <- dhalfcauchy(sigma.alpha, 25, log=TRUE)
    sigma.prior2 <- dhalfcauchy(sigma.beta, 25, log=TRUE)
    sigma.prior3 <- dhalfcauchy(sigma, 25, log=TRUE)

    log.prior <- mean.alpha.prior + mean.beta.prior + sum(alphas.prior) + sum(betas.prior) + sigma.prior1 + sigma.prior2 + sigma.prior3

    # predictions
    LogPreds <- alphas[Index] + Betas[Index]*log(Lengths)

    # log-Likelihood
    LLComp <- dnorm(log(Weights),LogPreds,sigma,log=TRUE)
    LL <- sum(LLComp)

    LP <- LL + log.prior
    Modelout <- list(LP=LP,Dev=-2*LL, Monitor= c(LP,sigma,sigma.alpha,sigma.beta), yhat=exp(rnorm(length(LogPreds), LogPreds, sigma)),parm=parm)
    show(proc.time() - a)
    return(Modelout)
  }

  Initial.Values <- c(log(0.01),log(0.05),log(3),log(0.1),rep(log(0.01),Npop),rep(log(3),Npop),log(1))

  andres <- Model(Initial.Values,Data = MyData)

  #print(Initial.Values)
  Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,Covar=NULL, Iterations=1e6,
                       Status=10000, Thinning=20000,Algorithm="HARM",Specs=list(alpha.star=0.23, B = NULL))
#   plot(Fit,BurnIn=100,MyData,PDF=FALSE,Parms=NULL)
#
#   Pred <- predict(Fit,Model,MyData)
#   Pred$y <- Weights
#   plot(Pred,Style="Fitted")
#   plot(Pred,Style="Residuals")
#
#   print(Fit)
  return(Fit)
}

GenerateData()
TheData <- read.table("G:\\Data.DAT",header=F)
Index <- TheData[,1]
Lengths <- TheData[,2]
Weights <- TheData[,3]

Fit <- Analyze()
Consort(Fit)
