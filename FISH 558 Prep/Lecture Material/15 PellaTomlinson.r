##############################################
#Plotting the Pella-Tomlinson model with 
#z=1.188 so that MSY occurs at 40% of B0, then
#reformulated so that it can handle B/Bmsy values
#Used in Branch, T.A., Hively, D.J., and Hilborn, R. 2013. 
#Is the ocean food provision index biased? Nature 495: E5-E6.
#This R code prepared for FISH 458
#Trevor A Branch, tbranch@uw.edu
#starting 2 May 2013
##############################################

pella.tomlinson.branch <- function() {
   z <- 1.188
   BdivBmsy <- seq(0,3.0,length.out=1000)
   FPI <- vector(length=length(BdivBmsy))
   
   for(i in 1:length(BdivBmsy)) {
      temp <- (1/z)^(1/(z-1))
      FPI[i] <- (z^(z/(z-1)))/(z-1) * (temp*BdivBmsy[i] - ((temp*BdivBmsy[i])^z))
      FPI[i] <- max(0, FPI[i])
   }
   plot(x=BdivBmsy, y=100*FPI, xaxs="i", yaxs="i", ylim=c(0,105), col="darkgreen",lwd=2.5, las=1,
        xlab="B/Bmsy", ylab="Food provision (%)", type="l")
}
 pdf("Figs\\PT model checking Halpern.pdf",width=8, height=5)
 
pella.tomlinson.branch()
 dev.off()