Ex3Soln <- function()
{
 DoSir(Nout=200,Model=F)
 DoSir(Nout=200,Model=T)

}

# =================================================================================

DoSir <- function(Nout=1000,Model)
{

 TheData <- ReadData()
 Yr1 <- TheData$CatchYr[1]
 Catch <- TheData$CatchVal
 SurveyEst <- TheData$SurveyEst
 SurveyCV <- TheData$SurveyCV
 SurveyYr <- TheData$SurveyYr
 Nyears <- length(Catch)

 Threshold <- exp(0)
 Cumu <- 0
 Ndone <- 0
 AboveK <- 0
 Vals <- matrix(0,ncol=5,nrow=Nout)
 AveLike <- 0
 Ntest <- 0
 while (Ndone < Nout)
  {
   r <- runif(1,0,0.15)
   K <- runif(1,20000,50000)
   Pop1965 <- runif(1,10000,15000)
   AddCV <- runif(1,0.1,0.2)
  # r <- 0.05
   #K <- 20000
   #Pop1965 <- 14000
   #AddCV <- 0.1
   #Model <- T
   Pop <- PopModel(Catch,r,K,Pop1965,Model)
   NegLogLike <- Likelihood(Pop,SurveyYr-Yr1+1,SurveyEst,SurveyCV,AddCV)
   TheLike <- exp(-1*NegLogLike-32.19)
   Cumu <- Cumu + TheLike
   AveLike <- AveLike + TheLike
   Ntest <- Ntest + 1
   while (Cumu > Threshold & Ndone < Nout)
    {
     cat("Saving",Ndone,"\n")
     Cumu <- Cumu - Threshold
     Ndone <- Ndone + 1
     Vals[Ndone,1] <- r
     Vals[Ndone,2] <- K
     Vals[Ndone,3] <- Pop1965
     Vals[Ndone,4] <- AddCV
     Vals[Ndone,5] <- Pop[Nyears]/K
     if (max(Pop/K) > 0.9) AboveK <- AboveK + 1/Nout
    }
  }
 cat(AboveK,AveLike/Ntest,Ntest,"\n")

 par(mfrow=c(3,2))
 hist(Vals[,1],main="",xlab="r")
 hist(Vals[,2],main="",xlab="K")
 hist(Vals[,3],main="",xlab="Population size 1965")
 hist(Vals[,4],main="",xlab="Addition CV")
 hist(Vals[,5],main="",xlab="Depletion")


}

# =================================================================================

Likelihood <- function(Pop,SurveyYr,SurveyEst,SurveyCV,AddCV)
{
 UseCV <- sqrt(SurveyCV^2+AddCV^2)
 Preds <- Pop[SurveyYr]
 Residuals <- log(UseCV)+0.5*(log(Preds)-log(SurveyEst))^2/UseCV^2
 LogLike <- sum(Residuals)
}

# =================================================================================
PopModel <- function(Catch,r,K,InitPop,ExponModel)
{
 Nyears <- length(Catch)
 Pop <- rep(0,length=Nyears)

 Pop[1] <- InitPop
 for (Iyear in 1:(Nyears-1))
  if (ExponModel == T)
   Pop[Iyear+1] <- max(0.01,Pop[Iyear]*(1+r)-Catch[Iyear])
  else
   Pop[Iyear+1] <- max(0.01,Pop[Iyear] + r*Pop[Iyear]*(1-Pop[Iyear]/K) - Catch[Iyear])
 return(Pop)

}

# =================================================================================
ReadData <- function()
{
 FileName <- "G:\\Ex3a.csv"
 TheData1 <- matrix(scan(FileName,skip=1,n=22*3,sep=','),ncol=3,byrow=T)
 FileName <- "G:\\Ex3b.csv"
 TheData2 <- matrix(scan(FileName,skip=1,n=38*2,sep=','),ncol=2,byrow=T)

 Outs <- NULL
 Outs$SurveyYr <- TheData1[,1]
 Outs$SurveyEst <- TheData1[,2]
 Outs$SurveyCV <- TheData1[,3]
 Outs$CatchYr <- TheData2[,1]
 Outs$CatchVal <- TheData2[,2]
 return(Outs)
}

# =================================================================================
Ex3Soln()
