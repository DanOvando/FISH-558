####FISH 458 lecture on spatial models
#Written by Trevor A. Branch starting 6 May 2012
#Completed 7 May 2012
#######################################################

###Model 1: one-D spatial model logistic growth with harvesting
###assume reflection when hits boundary
###Test to see what total yield is at different harvest rates, and with 
#different numbers of cells closed

one.d.logistic <- function(ncells=101, expl.rate=rep(0,ncells), ntime.steps=20, K=1000, r=0.2, m.rate=0.1,
                           start.N.vec=c(rep(0,50),K,rep(0,50))) {
  N.mat <- matrix(nrow=ntime.steps,ncol=ncells)
  
  N.mat[1,] <- start.N.vec
  yield <- vector(length=ntime.steps)
  yield[] <- 0
  for (year in 1:(ntime.steps-1)) {
    yield[year] <- 0
    for (i in 1:ncells) {
       yield[year] <- yield[year]+expl.rate[i]*N.mat[year,i]
    }
    #modelled as reflection, so no movement when hit bounds
    #left-most cell
    N.mat[year+1,1] <- N.mat[year,1] + r*N.mat[year,1]*(1-N.mat[year,1]/K) -
          expl.rate[1]*N.mat[year,1] +
          m.rate*((1-expl.rate[1+1])*N.mat[year,1+1]) -
          m.rate*((1-expl.rate[1])*N.mat[year,1])

    #right-most cell
    N.mat[year+1,ncells] <- N.mat[year,ncells] + r*N.mat[year,ncells]*(1-N.mat[year,ncells]/K) -
          expl.rate[ncells]*N.mat[year,ncells] +
          m.rate*((1-expl.rate[ncells-1])*N.mat[year,ncells-1]) -
          m.rate*((1-expl.rate[ncells])*N.mat[year,ncells])
    
    #all the cells in between
    for (i in 2:(ncells-1)) {
      N.mat[year+1,i] <- N.mat[year,i] + r*N.mat[year,i]*(1-N.mat[year,i]/K) -
          expl.rate[i]*N.mat[year,i] + 
          m.rate*((1-expl.rate[i-1])*N.mat[year,i-1]+(1-expl.rate[i+1])*N.mat[year,i+1]) -
          2*m.rate*((1-expl.rate[i])*N.mat[year,i])
      #print(N.mat[year+1,i])
    }
  }
  #now calculate final year yield
  yield[ntime.steps] <- 0
  for (i in 1:ncells) {
    yield[ntime.steps] <- yield[ntime.steps]+expl.rate[i]*N.mat[ntime.steps,i]
  }
  return(list(N.mat=N.mat, yield=yield))
}


#==============================calculations and plotting============================
####show the time series of yield by time for a few exploitation rates
K.start <- 1000
ncells <- 51
ntime.steps <- 100    #use 1000 for plotting, 100 for class demo
urate <- c(0.1,0.25)
nharvest <- length(urate)
yield.time.series <- matrix(nrow=ntime.steps,ncol=nharvest)
for (i in 1:nharvest) {
   x <- one.d.logistic(ncells=ncells, 
                       expl.rate=c(rep(urate[i],(ncells-1)/2-2),0,0,0,0,0,rep(urate[i],(ncells-1)/2-2)), 
                       ntime.steps=ntime.steps, 
                       K=K.start, r=0.2, m.rate=0.2,
                       start.N.vec=rep(K.start,ncells)) 
   yield.time.series[,i] <- x$yield
}

yield.time.series.noMPA <- matrix(nrow=ntime.steps,ncol=nharvest)
for (i in 1:nharvest) {
   xnoMPA <- one.d.logistic(ncells=ncells, 
                            expl.rate=rep(urate[i],ncells), ntime.steps=ntime.steps, 
                            K=K.start, r=0.2, m.rate=0.2,
                            start.N.vec=rep(K.start,ncells)) 
   yield.time.series.noMPA[,i] <- xnoMPA$yield
}

plot(x=1,y=1, xaxs="i",yaxs="i",type="n",lwd=2,col="blue",las=1,
     xlab="Year", ylab="Yield in each year", ylim=c(0,13000),cex.lab=1.3, xlim=c(0,ntime.steps))

lines(x=1:ntime.steps,y=yield.time.series[,1],type="l",lwd=2,col="blue",las=1)
lines(x=1:ntime.steps,y=yield.time.series[,2],type="l",lwd=2,col="blue",las=1)
lines(x=1:ntime.steps,y=yield.time.series.noMPA[,1],type="l",lwd=2,col="green",las=1)
lines(x=1:ntime.steps,y=yield.time.series.noMPA[,2],type="l",lwd=2,col="green",las=1)


####for demo use ntime.steps=100, for plotting use ntimes.steps=1000

###compare yield by harvest rates for 5 closed cells out of 51
K.start <- 1000
ncells <- 51
ntime.steps <- 100    #use 1000 for plotting
urate <- seq(0,0.2,0.01)
nharvest <- length(urate)
final.yield <- vector(length=nharvest)
for (i in 1:nharvest) {
  x <- one.d.logistic(ncells=ncells, 
                      expl.rate=c(rep(urate[i],(ncells-1)/2-2),0,0,0,0,0,rep(urate[i],(ncells-1)/2-2)), 
                      ntime.steps=ntime.steps, 
                      K=K.start, r=0.2, m.rate=0.2,
                      start.N.vec=rep(K.start,ncells)) 
  final.yield[i] <- x$yield[ntime.steps]
}

final.yield.noMPA <- vector(length=nharvest)
for (i in 1:nharvest) {
  xnoMPA <- one.d.logistic(ncells=ncells, 
                      expl.rate=rep(urate[i],ncells), ntime.steps=ntime.steps, 
                      K=K.start, r=0.2, m.rate=0.2,
                      start.N.vec=rep(K.start,ncells)) 
  final.yield.noMPA[i] <- xnoMPA$yield[ntime.steps]
}

plot(x=urate,y=final.yield, xaxs="i",yaxs="i",type="l",lwd=2,col="blue",las=1,
           xlab="Harvest rate", ylab="Equilibrium yield", ylim=c(0,2600),cex.lab=1.3)
par(new=T)
plot(x=urate,y=final.yield.noMPA, xaxs="i",yaxs="i",type="l",lwd=2,col="green",las=1,
     xlab="", ylab="",ylim=c(0,2600),axes=F)


###compare yield by harvest rates for 25 closed cells out of 51
K.start <- 1000
ncells <- 51
ntime.steps <- 1000
urate <- seq(0,0.9,0.01)
nharvest <- length(urate)
final.yield <- vector(length=nharvest)
for (i in 1:nharvest) {
  x <- one.d.logistic(ncells=ncells, 
                      expl.rate=c(rep(urate[i],13),rep(0,25),rep(urate[i],13)), ntime.steps=ntime.steps, 
                      K=K.start, r=0.2, m.rate=0.2,
                      start.N.vec=rep(K.start,ncells)) 
  final.yield[i] <- x$yield[ntime.steps]
}

final.yield.noMPA <- vector(length=nharvest)
for (i in 1:nharvest) {
  xnoMPA <- one.d.logistic(ncells=ncells, 
                           expl.rate=rep(urate[i],ncells), ntime.steps=ntime.steps, 
                           K=K.start, r=0.2, m.rate=0.2,
                           start.N.vec=rep(K.start,ncells)) 
  final.yield.noMPA[i] <- xnoMPA$yield[ntime.steps]
}

plot(x=urate,y=final.yield, xaxs="i",yaxs="i",type="l",lwd=2,col="blue",las=1,
     xlab="Harvest rate", ylab="Equilibrium yield", ylim=c(0,2600),cex.lab=1.3)
par(new=T)
plot(x=urate,y=final.yield.noMPA, xaxs="i",yaxs="i",type="l",lwd=2,col="green",las=1,
     xlab="", ylab="",ylim=c(0,2600),axes=F)


###compare yield by harvest rates for 25 closed cells out of 51
###Three scenarios that produce identical yield of around 1377 tons
###1. closed 25 areas, at u = 0.1
###2. all areas open at u = 0.032
###3. all areas open at u = 0.168

K.start <- 1000
ncells <- 51
ntime.steps <- 1000
urate <- 0.1
x <- one.d.logistic(ncells=ncells, 
                      expl.rate=c(rep(urate,13),rep(0,25),rep(urate,13)), ntime.steps=ntime.steps, 
                      K=K.start, r=0.2, m.rate=0.2,
                      start.N.vec=rep(K.start,ncells)) 
print(x$yield[ntime.steps])

urate<-0.032
xnoMPA1 <- one.d.logistic(ncells=ncells, 
                           expl.rate=rep(urate,ncells), ntime.steps=ntime.steps, 
                           K=K.start, r=0.2, m.rate=0.2,
                           start.N.vec=rep(K.start,ncells)) 
print(xnoMPA1$yield[ntime.steps])

urate<-0.168
xnoMPA2 <- one.d.logistic(ncells=ncells, 
                          expl.rate=rep(urate,ncells), ntime.steps=ntime.steps, 
                          K=K.start, r=0.2, m.rate=0.2,
                          start.N.vec=rep(K.start,ncells)) 
print(xnoMPA2$yield[ntime.steps])

plot(x=1:ncells,y=x$N.mat[ntime.steps,], ylim=c(0,1050), xlab="Cell number",cex.lab=1.3,ylab="Abundance",type="l",lwd=2,col="blue", las=1,xaxs="i",yaxs="i")
par(new=T)
plot(x=1:ncells,y=xnoMPA1$N.mat[ntime.steps,], ylim=c(0,1050), xlab="",cex.lab=1.3,ylab="",type="l",lwd=2,col="green", las=1,xaxs="i",yaxs="i",axes=F)
par(new=T)
plot(x=1:ncells,y=xnoMPA2$N.mat[ntime.steps,], ylim=c(0,1050), xlab="",cex.lab=1.3,ylab="",type="l",lwd=2,col="green", las=1,xaxs="i",yaxs="i",axes=F)
polygon(x=c(14,14,38,38),y=c(0,1050,1050,0),col="#77777733")
