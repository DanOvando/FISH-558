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
#this function returns a vector of harvest rates in each cell. 
#===================================================================
calc.harvest.vec <- function(u.out, ncells, MPA.width)  {

}
calc.harvest.vec(u.out=0.1, ncells=4, MPA.width=2)

#======PART 2====================================================
#Creates a simple MPA model with logistic growth
#r=rate of increase, K=carrying capacity, u.out=harvest rate outside MPA,
#ncells=number of cells for model, nsteps=number of time periods,
#MPA.width=number of cells in MPA, mrate=movement rate
#===================================================================
MPA.model <- function(r, K, u.out, ncells, nsteps, MPA.width, mrate) {
   u.vec <- calc.harvest.vec(u.out=u.out, ncells=ncells, MPA.width=MPA.width)
   pop <- vector(length=ncells)

   #set starting population in vector pop in each cell equal to K

   
   #vector left.cells storing the position of cells to the left

   
   #vector right.cells storing the position of cells to the right

   
   
   #plot the initial numbers
   plot(x=1:ncells, y=pop, xlab="Cell number", lwd=3,
        ylab="Population size", ylim=c(0, 1.05*max(pop)), type="l", yaxs="i", xaxs="i")
   

   #loop through the time steps
   for (i in 1:nsteps) {
      #Vector "leaving" of number leaving each cell is population size times diffusion rate

      
      #The number of immigrants is 1/2 those leaving cells to the
      #left and 1/2 those leaving cells to the right. This is a
      #complicated expression, take a minute to think it through!
      arriving <- 0.5*leaving[left.cells]+0.5*leaving[right.cells]
      
      #surplus production from the logistic model

      
      
      #catches = harvest rate in each cell times the population size

      
      
      #update the population numbers

      
      #plot the population in each cell
      #lines(x=1:ncells, y=pop, lwd=(nsteps-i+1)/nsteps*3)
   }
}
MPA.model(r=0.2, K=1000, u.out=0.4, ncells=21, MPA.width=5, nsteps=10, mrate=0.4)


#======PART 3====================================================
#Equilibrium sum of catches. 
#1. Copy function from Part 2.  
#2. Remove the plotting functions (so it runs faster)
#3. Return sum of catches after 1000 time steps.
#===================================================================
MPA.eqm.catch <- function(r, K, u.out, ncells, nsteps, MPA.width, mrate) {

   
   

}
MPA.eqm.catch(r=0.2, K=1000, u.out=0.4, ncells=21, MPA.width=5, nsteps=1000, mrate=0.4)

#======PART 4====================================================
#Plots contours of catches on a plot of MPA width (y axis) vs 
#harvest rate (x axis).
#r=rate of increase, K=carrying capacity, 
#ncells=number of cells for model, nsteps=number of time periods,
#mrate=movement rate
#===================================================================
contour.width.harvest <- function(r=0.2, K=1000, ncells=100, nsteps=1000, mrate=0.4) {
   #vector MPA.widths[] to loop over 
   
   
   #vector harvest[] to loop over from 0*r to 1*r

   
   
   #matrix catch.mat[] to store the equilibrium catches
   
   
    
   #loop through widths and loop through harvest rates
   #storing equilm summed catches from MPA.eqm.catch() function in the catch matrix

   
   
   
   
   #Plot a contour showing total catches by MPAwidth and harvest rate
   #xvalues are the harvest values, y values are MPA.widths, and 
   #z values are equilibrium catches stored in the matrix
   contour(x=harvest, y=MPA.widths, z=catch.mat,
         xlab="Harvest rate",  ylab="MPA width", xaxs="i",yaxs="i")
   
   #return the resulting values for catches
   return(catch.mat)
}
x <- contour.width.harvest(r=0.2, K=1000, ncells=100, nsteps=1000, mrate=0.2) 
x
