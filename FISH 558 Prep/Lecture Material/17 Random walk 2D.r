####Demonstration of random walk in 2 dimensions of space
####Written by Trevor A. Branch starting 7 May 2012
####tbranch@uw.edu

random.2D <- function(start.x = 0, start.y = 0, nsteps=50, step.size=0.5) {
    locs <- matrix(nrow=nsteps,ncol=3)
    locs[,1] <- 1:nsteps   #first column stores the years (time steps)
    locs[1,2:3] <- c(start.x,start.y)  #starting position in time step 1
    for (yr in 1:(nsteps-1)) {
       #pick two random numbers between -step and +step and add them to the x and y
       locs[yr+1,2:3] <- locs[yr,2:3]+runif(n=2,min=-step.size, max=step.size)
    }
    #find the square the encloses the results
    xxx <- max(abs(min(locs[,2])),max(locs[,2]))
    yyy <- max(abs(min(locs[,3])),max(locs[,3]))
    
    #create a blank plot with the correct x limits and y limits
    plot(type="n",x=1,axes=F,xlab="",ylab="",xlim=c(-xxx,xxx),
         ylim=c(-yyy,yyy))
    
    #use segments to plot lines that change in color over time
    segments(x0=locs[-1,2],x1=locs[-nsteps,2],y0=locs[-1,3],y1=locs[-nsteps,3],
       col=rev(gray((1:nsteps / nsteps)*0.8)),lwd=2)
    #axis(1)  #if you want coordinates uncomment this line
    #axis(2)
    abline(h=0,lty=2,col="grey50")
    abline(v=0,lty=2,col="grey50")
    
}
#pdf("Figs\\17 random walk2D.pdf", width=6,height=6)
par(mfrow=c(2,2),oma=c(0,0,0,0),mar=rep(0.5,4))
for (j in 1:4) {
  random.2D(start.x=0,start.y=0,nsteps=100, step.size=0.5)
}
#dev.off()
