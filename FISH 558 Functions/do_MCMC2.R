DoMCMC<-function(Xinit,DataUsed,Ndim,covar,Nsim=1000,Nburn=0,Nthin=1,prog_bar = T,run_time)
{

  Xcurr <- Xinit
  Fcurr <- -1*NegLogLike2(Xcurr = Xcurr,dat = DataUsed,run_time = run_time)
  Outs <- matrix(0,nrow=(Nsim-Nburn),ncol=(Ndim+1))
  Ipnt <- 0; Icnt <- 0

  if(prog_bar == T)
  {
    pb <- txtProgressBar(min = 1, max = Nsim, style = 3)
  }

  for (Isim in 1:Nsim)
  {
    NewX <- rmvnorm(1, mean=Xcurr$guess, sigma=covar)
    Xnext <- Xcurr
    Xnext$guess <- t(NewX)

    Fnext <- -1*NegLogLike2(Xcurr = Xnext,dat = DataUsed,run_time = run_time )
    if(prog_bar == T)
    {

      setTxtProgressBar(pb, Isim)

    }

    Rand1 <- log(runif(1,0,1))
    if (Fnext > Fcurr+Rand1)
    {Fcurr <- Fnext; Xcurr <- Xnext }
    if (Isim > Nburn & Isim %% Nthin == 0)
    {
      Icnt <- Icnt + 1
      Outs[Icnt,] <- c(Xcurr$guess,Fcurr)
    }
  }
  xx <- seq(1:Icnt)
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
  return(Outs[1:Icnt,])
}