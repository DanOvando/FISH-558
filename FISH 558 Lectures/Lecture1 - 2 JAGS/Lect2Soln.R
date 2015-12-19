library(runjags)
# setwd("C:\\courses\\FISH 558_15\\Lectures\\Bayes Workshop FHL\\Exercises\\")

ResultsFile <- "Lect2Res.dat"

model <- "model {
 for(i in 1 : N){
  Count[i] ~ dpois(true.y[Plot[i]]*Effort[i])
 }
 for(j in 1 : Nplot) {
   true.y[j] ~ dlnorm(mean,precision)
 }
 mean ~ dunif(-1000,1000)
 meanN <- exp(mean)
 sigma <- 1.0/sqrt(precision)
 precision ~ dunif(0,1000.0)
 }"

Cases <- c(1,0,1,0)

# Read in the data
# ================
TheD <- read.table("Lect2.dat",head=T)
Plot <- as.double(TheD$Plot)
Effort <- as.double(TheD$Effort)
Count <-  as.double(TheD$Count)
N <-length(Plot)
Nplot <-length(unique(Plot))

if (Cases[1] ==1)
 {
  data <- list(Plot=Plot, Effort=Effort, Count=Count,N=N,Nplot=Nplot)
  inits1 <- list(mean=log(20), precision=0.1,true.y=rep(20,Nplot),
                 .RNG.name="base::Super-Duper", .RNG.seed=1)
  inits2 <- list(mean=log(10), precision=0.2,true.y=rep(10,Nplot),
                 .RNG.name="base::Super-Duper", .RNG.seed=1)
  inits3 <- list(mean=log(30), precision=0.3,true.y=rep(30,Nplot),
                 .RNG.name="base::Super-Duper", .RNG.seed=1)
  inits <- list(inits1,inits2,inits3)

  browser()

  results <- run.jags(model=model, monitor=c("mean", "meanN","sigma","precision","true.y"),
                    data=data,n.chains=3, method="rjags", inits=inits,
                    plots=T,silent.jag=F,
                    sample=2000,adapt=2000,burnin=10000,thin=10)

  print(results)
  write(t(results$mcmc[[1]]),ResultsFile,ncolumn=14)
  write(t(results$mcmc[[2]]),ResultsFile,ncolumn=14,append=T)
  write(t(results$mcmc[[3]]),ResultsFile,ncolumn=14,append=T)

  TheRes <- read.table(ResultsFile)

  par(mfrow=c(2,2))
  d <- density(TheRes[,2])
  plot(d,xlab="Mean",main="")
  d <- density(TheRes[,3])
  plot(d,xlab="Sigma",main="")

}

# Posterior distribution for the mean density by plot
if (Cases[2]==1)
 {
  TheRes <- read.table(ResultsFile)
  par(mfrow=c(4,3),mar=c(5,4,1,1),omi=c(0.5,0.5,0.5,0.5))
  for (Iplot in 1:Nplot)
   {
    hist(TheRes[,5+Iplot],xlab=paste("Plot", Iplot),main="")
   }
 }

# posterior predictive distribution by data point

pdf(file = 'stupidjags.pdf')
if (Cases[3]==1)
 {
  TheRes <- read.table(ResultsFile)
  Nsim <- length(TheRes[,1])
  print(head(TheRes))
  par(mfrow=c(5,5),mar=c(5,4,1,1),omi=c(0.5,0.5,0.5,0.5))
  TheRes <- read.table(ResultsFile)
  for (Idata in 1:N)
   {
    Iplot <- Plot[Idata]+4
    Eff <- Effort[Idata]
    Obs <- Count[Idata]
    Vec <- rep(NA,Nsim)
    for (Isim in 1:Nsim)
     Vec[Isim] <- rpois(1,TheRes[Isim,Iplot]*Eff)
    hist(Vec,xlab=paste("Data point",Idata),main="")
    abline(v=Obs,lwd=4)
   }
}
dev.off()

# posterior predictive distribution by data point
if (Cases[4]==1)
 {
  TheRes <- read.table(ResultsFile)
  Nsim <- length(TheRes[,1])
  Nsim <- length(TheRes[,1])
  Effs <- c(5,10,15)
  resu <- matrix(NA,nrow=Nsim,ncol=3)
  for (Isim in 1:Nsim)
   {
    Pdens <- rlnorm(1,TheRes[Isim,2],TheRes[Isim,5])
    for (Ieff in 1:3)
     resu[Isim,Ieff] <- rpois(1,Pdens*Effs[Ieff])
   }
  par(mfrow=c(2,2))
   for (Ieff in 1:3)
    hist(resu[,Ieff],xlab=paste("Effort =",Effs[Ieff]),main="")
}
