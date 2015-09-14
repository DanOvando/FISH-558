#================================================================
#Lab 5 for FISH458: PART 4 adding predator functional response
#Building a predator-prey model in R
#Using the sliders in Rstudio to explore the effect of changing 
#model parameters.
#Written by Trevor A. Branch, tbranch@uw.edu
#Started on 29 April 2013
#================================================================
predator.prey <- function(wild0, wildK, lion0, killrate, assimil, nyears, 
                          wildr=0.2, lionsurv=0.8,  
                          arrow.length=0.0, 
                          aprime=1460, pc=0.001, h=0.01369863, Tt=1, TotalArea=90000) {
   #define variables to store the numbers of wildebeest and lions in each year
   N.wild <- vector(length=nyears)
   N.lions <- vector(length=nyears)
   
   #set the numbers of wildebeest and lions in the first year
   N.wild[1] <- wild0
   N.lions[1] <- lion0
   
   #calculate the number of wildebeest and lions in each year
   for (i in 1:(nyears-1)) {
      Nt <- N.wild[i]/TotalArea
      X <- Tt*aprime*pc*Nt/(1+h*aprime*pc*Nt)
      deaths <- N.lions[i]*X
      N.wild[i+1] <- N.wild[i]+wildr*N.wild[i]*(1-N.wild[i]/wildK)-deaths
      N.lions[i+1] <- lionsurv*N.lions[i]+assimil*deaths
   }
   
   #make some space below and to the left for the axes, create a 2x1 plot area
   par(mfrow=c(2,1), oma=c(0,0,0,0), mar=c(5,5,1,5))
   
   #plot numbers by year
   plot(x=1:nyears, y=N.wild/1000, type="l", lwd=2, col="blue",
        ylab="", xlab="Years", cex.lab=1.2, las=1)
   par(new=T)  #plot on the same plot
   plot(x=1:nyears, y=N.lions/1000, type="l", lwd=2, col="red", yaxt="n", xaxt="n", 
        ylab="", xlab="Years", cex.lab=1.2)   #suppress axis
   axis(side=4, las=1)  #put another axis on the right for lions
   mtext(side=2, "Wildebeest (in thousands)", line=3, cex=1.1, col="blue")
   mtext(side=4, "Lions (in thousands)", line=2.5, cex=1.1, col="red")
   
   #phase plane, first create an empty plot with axes 5% bigger than the data and labels for the axes
   plot(x=1, y=1, type="n", xaxs="i", yaxs="i", las=1, 
        xlim=c(0,1.05*max(N.wild))/1000, ylim=c(0,1.05*max(N.lions))/1000, 
        ylab="Lions (in thousands)", xlab="Wildebeest (in thousands)", cex.lab=1.2)
   #plot a series of arrows getting darker with the years
   #to show arrows, make arrow.length=0.07
   arrows(x0=N.wild[-nyears]/1000, x1=N.wild[-1]/1000, y0=N.lions[-nyears]/1000, y1=N.lions[-1]/1000, lwd=2, 
               length=arrow.length, col=gray(0.9*(nyears:1)/nyears))
}
predator.prey(wild0=1000000, lion0=8000, wildr=0.2, wildK=1500000, lionsurv=0.8, 
              killrate=0.000016, assimil=0.05, nyears=300, arrow.length=0.0)

library(manipulate)
manipulate(predator.prey(wild0, wildK, lion0, killrate, assimil, nyears), 
           wild0=slider(100000, 2000000, initial=1000000), 
           wildK=slider(1000000, 6000000, initial=1500000), 
           lion0=slider(100, 10000, initial=8000),
           killrate=slider(0.000005, 0.000035, initial=0.000016),
           assimil=slider(0.01, 0.2, initial=0.05),
           nyears=slider(10,1000, initial=300))