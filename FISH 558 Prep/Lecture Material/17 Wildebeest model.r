################################################################################
###Wildebeest foraging model designed by Ray Hilborn in Visual Basic in Excel
##Translated to R by Trevor A. Branch tbranch@gmail.com starting 8 May 2012
##The animation output requires installation of ImageMagick, www.imagemagick.org/script/index.php, followed by machine restart
##assumes grass is not consumed by wildebeest
################################################################################
run.model <- function(Xbar=5.1, Ybar=5.1, Xsigma=0.11,Ysigma=0.11,nX=10,nY=10,dBar=1,
                       nWild=50,nYear=40, strong=T, movie.file="wildebeest40strong.gif", 
                      outdir ="C:\\Users\\Trevor Branch\\Documents\\2013 FISH458 Quantitative\\1 Lectures\\Rpics") {
  #set up animation
  require(animation)
  oopts = ani.options(interval = 0.03, ani.dev="png", outdir = outdir, width=1400, height=600)   #interval is wait in seconds between frames in video
  par(bg = "white")        #ensure the background color is white
  ani.record(reset = TRUE) #clear history before recording
  
  
  ####Set up all vectors and matrices
  xList <- vector(length=9)
  yList <- vector(length=9)
  Grass <- matrix(nrow=nX,ncol=nY)   #grass quality in each location of the matrix
  Wild <- matrix(nrow=nX,ncol=nY)    #number of wildebeest in each location of the matrix
  xLoc <- vector(length=nWild)       #x location of each individual wildebeest
  yLoc <- vector(length=nWild)       #y location of each wildebeest
  Grass <- matrix(nrow=nX,ncol=nY)
  
  ##initialize grass quality in each square
  for (i in 1:nX) {
    for (j in 1:nY) {
      if (strong==F) {
        Grass[i,j] = sqrt(dBar / sqrt(((i-Xbar)^2*Xsigma + (j-Ybar)^2*Ysigma)) ^ 2)  #initial grass abundance in each cell
      }
      else {
        Grass[i,j] = dBar / sqrt(((i-Xbar)^2*Xsigma + (j-Ybar)^2*Ysigma)) ^ 2  #initial grass abundance in each cell
      }
    }
  }
  
  ##initialize the location of wildebeest
  for (i in 1:nWild) {  #all wildebeest start in square 1
    xLoc[i] <- 1
    yLoc[i] <- 1
  }
  
  ##to convert from vector to 3x3 array
  xList  <- c(-1,  0,  1,             -1,  0,  1,             -1,  0,  1)
  yList  <- c(-1, -1, -1,              0,  0,  0,              1,  1,  1)
  
  ###Plot the distribution of grass and of wildebeest
  par(oma=rep(0,4),mar=c(3,3,1,1))
  image(x=1:nX,y=1:nX,z=log(Grass),col=rev(heat.colors(100)),axes=F,xlab="",ylab="")
  par(new=T)
  plot(x=jitter(xLoc,amount=0.4),y=jitter(yLoc,amount=0.4),type="p",col="black",pch=19,
       xlab="",ylab="",cex=0.9,xlim=c(0.5,nX+0.5),ylim=c(0.5,nY+0.5),xaxs="i",yaxs="i")
  box(col="gray70")
  ani.record()           #record the current frame
  
  #now loop over a bunch of time steps
  for (yr in 1:nYear) {
    #initialize number of wildebeest in each cell to 0
    for (i in 1:nX) {
      for (j in 1:nY) {
        Wild[i,j] <- 0      
      }
    }
    
    #now loop over individual wildebeest
    for (k in 1:nWild) {
      #Get the quality of cells surrounding x and y
      quality <- qual.func(iX=xLoc[k], iY=yLoc[k], xList=xList, yList=yList, Grass=Grass, nX=nX, nY=nY)  
      xRand <- runif(n=1)
      foundit <- F   #this is a flag
      for (i in 1:9) {
        if (xRand <= quality[i] & foundit==F) {   #remember quality is cumulative distribution
          xLoc[k] <- xLoc[k]+xList[i]
          yLoc[k] <- yLoc[k]+yList[i]
          Wild[xLoc[k],yLoc[k]] <- Wild[xLoc[k],yLoc[k]] + 1
          foundit <- T
        }  #end if 
      }  #end for i
    }  #end for k
  
    ###Plot the distribution of grass and of wildebeest
    image(x=1:nX,y=1:nX,z=log(Grass),col=rev(heat.colors(100)),axes=F,xlab="",ylab="")
    par(new=T)
    plot(x=jitter(xLoc,amount=0.4),y=jitter(yLoc,amount=0.4),type="p",col="black",pch=19,
         xlab="",ylab="",cex=0.9,xlim=c(0.5,nX+0.5),ylim=c(0.5,nY+0.5),xaxs="i",yaxs="i")
    box(col="gray70")
    ani.record()           #record the current frame
  }
  for (i in 1:7) {  ani.record()  }           #add a natural break record the current frame
  saveGIF(ani.replay(), img.name = "oneDmovie", movie.name=movie.file, convert="convert",
          ani.width=400, ani.height=400, interval = 0.5)    #ImageMagick
  
}

###function that computes the quality of surrounding cells
qual.func <- function(iX,iY,xList,yList, Grass, nX, nY) {
  tempQual <- vector(length=9)
  finalQual <- vector(length=9)
  
  sumQual <- 0
  for (i in 1:9) {
    if (((iX +xList[i]) %in% 1:nX) & ((iY +yList[i]) %in% 1:nY)) {
      tempQual[i] <- Grass[iX +xList[i], iY+yList[i]]
    }
    else {
      tempQual[i] <- 0  #outside the Serengeti
    }
    sumQual <-sumQual+tempQual[i]
  }
  hQual <- cumsum(tempQual)/sumQual
  return(hQual)
}
run.model(strong=T, movie.file="wildebeest40strong.gif")
run.model(strong=F, movie.file="wildebeest40nostrong.gif")
