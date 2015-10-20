library(mvtnorm)
library(stats4)
library(MASS)
library(coda)


# =================================================================================================================

DoMCMC<-function(Xinit,DataUsed,Ndim,covar,Nsim=1000,Nburn=0,Nthin=1)
{


  Xcurr <- Xinit
  Fcurr <- -1*NegLogLike2(Xcurr = Xcurr,dat = DataUsed)
  Outs <- matrix(0,nrow=(Nsim-Nburn),ncol=(Ndim+1))
  Ipnt <- 0; Icnt <- 0
  for (Isim in 1:Nsim)
  {
    NewX <- rmvnorm(1, mean=Xcurr$guess, sigma=covar)

    Xnext <- Xcurr
    Xnext$guess <- t(NewX)
    Fnext <- -1*NegLogLike2(Xcurr = Xnext,dat = DataUsed)
    Rand1 <- log(runif(1,0,1))
    if (Fnext > Fcurr+Rand1)
    {Fcurr <- Fnext; Xcurr <- Xnext }
    if (Isim %% Nthin == 0)
    {
      Ipnt <- Ipnt + 1
      if (Ipnt > Nburn) {
        Icnt <- Icnt + 1; Outs[Icnt,] <- c(Xcurr$guess,Fcurr); }
    }
  }
  xx <- seq(1,Icnt)

  pdf(file = 'MCMC Diag Plot.pdf')
  par(mfrow=c(4,4),mar=c(3,4,2,1))
  for (II in 1:(Ndim+1))
  {
    yy <- Outs[,II][1:Icnt]
    if (II <= Ndim)
      lab1 <- paste("Parameter ",II)
    else
      lab1 <- "Posterior density"
    if (II <= Ndim) yy <- exp(yy)
    plot(xx,yy,xlab="Cycle number",ylab=lab1,type='b',pch=16,cex=0.02)
  }
  par(mfrow=c(4,4),mar=c(3,4,2,1))
  for (II in 1:(Ndim+1))
  {
    yy <- Outs[,II][1:Icnt]
    if (II <= Ndim)
      lab1 <- paste("Parameter ",II)
    else
      lab1 <- "Posterior density"
    if (II <= Ndim) yy <- exp(yy)
    hist(yy,ylab=lab1,main="")
  }
  dev.off()
  return(Outs[1:Icnt,])
}

# ================================================================================
