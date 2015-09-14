#====Lab 7 MPAs========================================================
#Develops a working one-dimensional MPA model with wraparound (Pacman) 
#boundaries, and diffusion movement of fish to adjacent squares. 
#Developed by Trevor A. Branch, tbranch@uw.edu
#Based on previous code by Cole Monnahan (2012)
#Current version started 13 May 2013
#======================================================================

#======PART 1=============================================
#given the harvest rate outside the MPA, u.out, the total number
#of cells ncells, and the width of the MPA, MPA.width
#returns a vector of harvest rates in each cell. 
#===================================================================
calc.harvest.vec <- function(u.out, ncells, MPA.width)  {
   u.vec <- vector(length=ncells)
   u.vec[] <- u.out       #set every cell to u.out, equivalent to u.vec[1:ncells] <- u.out   
   if (MPA.width > 0) {   #no need to do this if there is no MPA! 
      MPA.begin <- round((ncells-MPA.width)/2)+1  #start cell of MPA
      MPA.end <- MPA.begin + MPA.width -1         #end cell of MPA
      u.vec[MPA.begin:MPA.end] <- 0        
   }
   return(u.vec)
}
calc.harvest.vec(u.out=0.1, ncells=4, MPA.width=1)

#======PART 2====================================================
#Creates a simple MPA model with logistic growth
#r=rate of increase, K=carrying capacity, 
#ncells=number of cells for model, nsteps=number of time periods,
#MPA.width=number of cells in MPA, mrate=movement rate
#===================================================================
MPA.model <- function(r, K, u.out, ncells, nsteps, MPA.width, mrate) {
   u.vec <- calc.harvest.vec(u.out=u.out, ncells=ncells, MPA.width=MPA.width)
   pop <- vector(length=ncells)
   pop[] <- K   #start population size at K in each cell
   
   #the position of cells to the left
   left.cells <- c(ncells,1:(ncells-1))
   #the position of cells to the right
   right.cells <- c(2:ncells,1)
   
   #plot the initial numbers
   plot(x=1:ncells, y=pop, xlab="Cell number", lwd=3,
        ylab="Population size", ylim=c(0, 1.05*max(pop)), type="l", yaxs="i", xaxs="i")
   
   for (i in 1:nsteps) {
      #Number leaving each cell is population size times diffusion rate
      leaving <- mrate*pop
      #The number of immigrants is 1/2 those leaving cells to the
      #left and 1/2 those leaving cells to the right. This is a
      #complicated expression, take a minute to think it through!
      arriving <- 0.5*leaving[left.cells]+0.5*leaving[right.cells]
      
      #surplus production from the logistic model
      surplus.prod <- r*pop*(1-pop/K)
      
      #catches = harvest rate in each cell times the population size
      catches <- u.vec*pop
      
      #update the population numbers
      pop <- pop + surplus.prod - leaving + arriving - catches
      
      lines(x=1:ncells, y=pop, lwd=(nsteps-i+1)/nsteps*3)
   }
}

MPA.model(r=0.2, K=1000, u.out=0.4, ncells=21, MPA.width=5, nsteps=10, mrate=0.4)

#======PART 3====================================================
#Equilibrium catches
#Remove the plotting functions (so it runs faster), return sum of 
#catches after 1000 time steps.
#Creates a simple MPA model with logistic growth
#r=rate of increase, K=carrying capacity, 
#ncells=number of cells for model, nsteps=number of time periods,
#MPA.width=number of cells in MPA, mrate=movement rate
#===================================================================
MPA.eqm.catch <- function(r, K, u.out, ncells, nsteps, MPA.width, mrate) {
   u.vec <- calc.harvest.vec(u.out=u.out, ncells=ncells, MPA.width=MPA.width)
   pop <- vector(length=ncells)
   pop[] <- K   #start population size at K in each cell
   
   #the position of cells to the left
   left.cells <- c(ncells,1:(ncells-1))
   #the position of cells to the right
   right.cells <- c(2:ncells,1)
   
   for (i in 1:nsteps) {
      #Number leaving each cell is population size times diffusion rate
      leaving <- mrate*pop
      #The number of immigrants is 1/2 those leaving cells to the
      #left and 1/2 those leaving cells to the right. This is a
      #complicated expression, take a minute to think it through!
      arriving <- 0.5*leaving[left.cells]+0.5*leaving[right.cells]
      
      #surplus production from the logistic model
      surplus.prod <- r*pop*(1-pop/K)
      
      #catches = harvest rate in each cell times the population size
      catches <- u.vec*pop
      
      #update the population numbers
      pop <- pop + surplus.prod - leaving + arriving - catches
      
   }
   return(sum(catches))
}
MPA.eqm.catch(r=0.2, K=1000, u.out=0.4, ncells=21, MPA.width=5, nsteps=1000, mrate=0.4)

#======PART 4====================================================
#Contours of catches on a plot of MPA width (y axis) vs 
#harvest rate (x axis).
#r=rate of increase, K=carrying capacity, 
#ncells=number of cells for model, nsteps=number of time periods,
#MPA.width=number of cells in MPA, mrate=movement rate
#===================================================================
contour.width.harvest <- function(r=0.2, K=1000, ncells=100, nsteps=1000, mrate=0.4) {
    MPA.widths <- seq(0,ncells,by=1)  
    harvest <- seq(0,1.1,length.out=100)*r
    
    #matrix to store the equilibrium catches
    catch.mat <- matrix(NA, nrow=length(harvest), ncol=length(MPA.widths))
    
    #loop through widths and harvest rates, storing equilm catches
    for (i in 1:length(harvest)) {
       for (j in 1:length(MPA.widths)) {
           catch.mat[i,j] <- MPA.eqm.catch(r=r, K=K, u.out=harvest[i], 
                                           ncells=ncells, MPA.width=MPA.widths[j], 
                                           nsteps=nsteps, mrate=mrate)    
       }
    }
    #Plot a contour showing total catches by MPAwidth and harvest rate
    contour(x=harvest, y=MPA.widths, z=catch.mat,
            xlab="Harvest rate",  ylab="MPA width", xaxs="i",yaxs="i")
    return(catch.mat)
}
x <- contour.width.harvest(r=0.2, K=1000, ncells=100, nsteps=1000, mrate=0.2) 
x
