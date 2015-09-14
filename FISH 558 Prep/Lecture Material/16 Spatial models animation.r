####FISH 458 lecture on spatial models
#Written by Trevor A. Branch starting 6 May 2012
#Completed 7 May 2012
#######################################################

###Model 1: one-D spatial model logistic growth with harvesting
###assume reflection when hits boundary

one.d.logistic <- function(ncells=101, expl.rate=rep(0,ncells), ntime.steps=20, K=1000, r=0.2, m.rate=0.1,
                           start.N.vec=c(rep(0,50),K,rep(0,50))) {
  N.mat <- matrix(nrow=ntime.steps,ncol=ncells)
  N.mat[1,] <- start.N.vec
  for (year in 1:(ntime.steps-1)) {
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
  
  #plot results
  gray.pal <- (1:ntime.steps)/ntime.steps*0.8    #returns values between 0.2 (light gray) and 1 (black)
  #print(gray.pal)
  plot(x=1:ncells,y=N.mat[1,], type="l",xaxs="i",yaxs="i",las=1, col=gray(gray.pal[1]), 
          ylim=c(0,K*1.02), xlab="", ylab="", lwd=2)
  for (year in 2:ntime.steps) {
     par(new=T)
     lines(x=1:ncells,y=N.mat[year,], col=gray(gray.pal[year]), ylim=c(0,K*1.02), lwd=2)
  }
  return(N.mat)
}

K.start <- 1000
ncells <- 21
par(oma=c(0,0,0,0), mar=c(5,5,1,1))
x <- one.d.logistic(ncells=ncells, expl.rate=rep(0,ncells), ntime.steps=10, K=K.start, r=0.2, m.rate=0.1,
               start.N.vec=c(rep(0,9),K.start/100,K.start,K.start/100,rep(0,9)))  


###Set up animation; all in middle cell, no harvest
#requires installation of ImageMagick, www.imagemagick.org/script/index.php, followed by machine restart
require(animation)

par(bg = "white")        #ensure the background color is white
ani.record(reset = TRUE) #clear history before recording
outdir <- "C:\\Users\\Trevor Branch\\Documents\\2013 FISH458 Quantitative\\1 Lectures\\Rpics"
oopts = ani.options(interval = 0.03, ani.dev="png", outdir = outdir, width=1400, height=600)   #interval is wait in seconds between frames in video

for (year in 2:50) {
  x <- one.d.logistic(ncells=ncells, expl.rate=rep(0,ncells), ntime.steps=year, K=K.start, r=0.2, m.rate=0.1,
                      start.N.vec=c(rep(0,9),K.start/100,K.start,K.start/100,rep(0,9))) 
  ani.record()           #record the current frame
}
saveGIF(ani.replay(), img.name = "oneDmovie", movie.name="Diffusion.gif", convert="convert",
           ani.width=1400, ani.height=600, interval = 0.2)    #ImageMagick

#####all at K, add harvest but MPA in middle cells
require(animation)
K.start <- 1000
ncells <- 21
outdir <- "C:\\temp"
oopts = ani.options(interval = 0.03, ani.dev="png", outdir = outdir, width=1400, height=600)   #interval is wait in seconds between frames in video

par(bg = "white")        #ensure the background color is white
ani.record(reset = TRUE) #clear history before recording

for (year in 2:50) {
  x <- one.d.logistic(ncells=ncells, expl.rate=c(rep(0.2,8),0,0,0,0,0,rep(0.2,8)), ntime.steps=year, K=K.start, r=0.2, m.rate=0.1,
                      start.N.vec=rep(K.start,ncells)) 
  ani.record()           #record the current frame
}
saveGIF(ani.replay(), img.name = "MPA", movie.name="MPAout.gif", convert="convert", interval=0.2,
        ani.width=1400, ani.height=600)    #ImageMagick

for(i in 1:100) {dev.off()}

###for(i in 1:100) {dev.off()}
