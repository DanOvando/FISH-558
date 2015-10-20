rm(list  = ls())
library(LaplacesDemon)

# setwd("G:\\")
TheData <- matrix(scan("Hake13.Dat"),ncol=3,byrow=3)
print(TheData)
Years <- TheData[,1]
Catch <- TheData[,2]
CPUE <- TheData[,3]

mon.names <- c("LL","sigma","rpar","Kpar")
parm.names <- as.parm.names(list(K=0,r=0))
PGF <- function(Data) return(c(rnormv(1,9,10),rnormv(1,-1,10)))
print(PGF)
print(parm.names)

MyData <- list(Nyear=length(Years), PGF=PGF, CPUE=CPUE, Years=Years, Catch=Catch,
     mon.names=mon.names, parm.names=parm.names,y=CPUE)

Model <- function(parm, Data)
 {

  K <- exp(parm[1])
  r <- exp(parm[2])

  Nyear <- Data$Nyear
  Catch <- Data$Catch
  CPUE <- Data$CPUE


  Bio <- rep(NA,length=Data$Nyear+1)

  Bio[1] <- K
  for (Year in 1:Nyear)
   {
    Bio[Year+1] <- Bio[Year] + r*Bio[Year]*(1-Bio[Year]/K) - Catch[Year]
    if (Bio[Year+1] < 0.01) Bio[Year+1] <- 0.01
   }

  qest <- 0; N <- 0
  for (Year in 1:Nyear)
   { qest <- qest + log(CPUE[Year]/Bio[Year]); N <- N + 1 }
  qest <- exp(qest / N)

  CPUEHat <- rep(NA,length=Data$Nyear)
  SS <- 0
  for (Year in 1:Nyear)
   {
    CPUEHat[Year] <-qest*Bio[Year]
    SS <- SS + (log(CPUE[Year]) - log(CPUEHat[Year]))^2
   }
  Sigma <- sqrt(SS/N)
  LL <- -1*(N*log(Sigma) + N/2.0)
  rpar <- r
  Kpar <- K
  #cat(r,K,Sigma,LL,"\n")

  Modelout <- list(LP=LL, Dev=-2*LL, Monitor=c(LL,Sigma,rpar,Kpar),yhat=rnorm(length(CPUEHat), CPUEHat, Sigma), parm=parm)
  return(Modelout)
 }

set.seed(666)

############################  Initial Values  #############################
Initial.Values <- c(log(4000),-1)

Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
    Covar=NULL, Iterations=1000000, Status=10000, Thinning=200,
    Algorithm="HARM")
print(Fit)
Consort(Fit)

plot(Fit,BurnIn=100,MyData,PDF=FALSE,Parms=NULL)
Pred <- predict(Fit,Model,MyData)
Pred$y <- CPUE
summary(Pred,Discrep="Chi-Square")

plot(Pred,Style="Density",Rows=1:9)
par(mfrow=c(2,2))
plot(Pred,Style="Fitted")
plot(Pred,Style="Residuals")





