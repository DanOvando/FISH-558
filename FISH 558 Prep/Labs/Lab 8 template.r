#======LAB 8============================================================
#Univ Washington FISH458
#Lab on dispersal kernel, matrix multiplication, and MPAs
#With added fleet dynamics
#Trevor A. Branch 21 May 2013
#Modified from code written by Cole Monnahan May 2012. 
#=======================================================================


#=======PART 1: dispersal kernel========================================
DispersalKernel <- function(sigma,ncell) {
   ## Given sigma and ncell, this function calculates a dispersion
   ## kernel matrix based on Normal density as a function of the
   ## distance of a row (i) to a column (j). It returns an ncell x
   ## ncell matrix with normalized rows.
   
   ## Loop through each cell (i,j) and calculate the probability of
   ## moving from i to j.
   dk <- matrix(nrow=ncell,ncol=ncell)
   for (i in 1:ncell)  {               # loop over cells:  from
      for (j in 1:ncell)  {           #loop over cells: to
         ## Get distance between cells either directly, or by
         ## wrapping around the ends
         
         ## Calculate the distance between two sites, which is the
         ## absolute value of the difference in indices.
         regdistance <- abs(i-j)
         ## wrap around distance when i is small and j is high
         wrapdistance1 <- abs(i+ncell-j)
         ## wrap around distance when j is small and i is high
         wrapdistance2 <- abs(i-ncell-j)
         distance <- min(regdistance,wrapdistance1,wrapdistance2)
         
         d <- dnorm(x=distance, mean=0, sd=sigma)
         dk[i,j] <- d             # store the density in the matrix
      }                            # end of columns
      ## Normalize the row so the numbers add up to 1.
      dk[i,] <- dk[i,]/sum(dk[i,])
   }                                   # end of rows
   return(dk)
}
#explore different values of sigma and ncell


#========PART 2: matrix multiplication===================================
#Practice understanding how matrix multiplication works in R
#========================================================================

#**create a 2-D matrix M using DispersalKernel(). This will be 
#the movement matrix. Use sigma



#**create abundance matrix N with ncell number of rows and nyear number of columns




#**set the abundance cell 3 in the first year to 1000



#**multiply abundance in the first year by the dispersal matrix



#**loop through the years and in each year store the result in the N matrix



#plot the abundance over time in each cell
plot(0,xlim=c(0,nyear+1), ylim=c(0,max(N)), type="n", xlab="Time", ylab="Abundance")
for(i in 1:ncol(N)) {
   lines(x=1:nyear, y=N[,i], col=rainbow(ncell)[i],lwd=2)
}

#shortcut to equilibrium biomass
M2 <- M
for (i in 1:nyear) {
   M2 <- M2 %*% M
}
M2
N[1,] %*% M2


#================PART 3===========================================
#updating Lab 7 MPA model with the new matrix method
#=================================================================

#harvest rate in each cell from last week
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


#MPA model with logistic growth
#r=rate of increase, K=carrying capacity, 
#ncells=number of cells for model, nsteps=number of time periods,
#MPA.width=number of cells in MPA, mrate=movement rate
MPA.model.lab8 <- function(r, K, u.out, ncells, nyear, MPA.width, sigma) {
   u.vec <- calc.harvest.vec(u.out=u.out, ncells=ncells, MPA.width=MPA.width)
   
   #**create a matrix N with nyear rows and ncells columns

   
   
   #**set the abundance in the first year to K in each cell

   
   
   #**calculate the M matrix here (only need do this once!)

   
   
   #plot the initial numbers
   plot(x=1:ncells, y=N[1,], xlab="Cell number", lwd=3,
        ylab="Population size", ylim=c(0, 1.05*max(N)), type="l", yaxs="i", xaxs="i")
   
   for (i in 1:(nyear-1)) {
      #** calculate surplus production from the logistic model

      
      
      #** calculate catch numbers from harvest rate

      
      
      #**calculate pre-moving abundance in each cell (store in temp vector)

      
      
      #** use matrix multiplication to disperse premove abundance to cells next year

      
      #plot the abundance in the next year
      lines(x=1:ncells, y=N[i+1,], lwd=(nyear-i+1)/nyear*3)
   }
}

MPA.model.lab8(r=0.2, K=1000, u.out=0.4, ncells=21, MPA.width=5, nyear=10, sigma=0.5)


#================PART 4===========================================
#Adding fleet dynamics
#=================================================================

#inside or outside the fished area?
IsFished <- function(ncells, MPA.width)  {
   is.fished <- vector(length=ncells)
   is.fished[] <- 1  
   if (MPA.width > 0) {   #no need to do this if there is no MPA! 
      MPA.begin <- round((ncells-MPA.width)/2)+1  #start cell of MPA
      MPA.end <- MPA.begin + MPA.width -1         #end cell of MPA
      is.fished[MPA.begin:MPA.end] <- 0        
   }
   return(is.fished)
}
temp <- IsFished(ncells=10, MPA.width=3)


#This function allocates boats proportional to abundance
#outside the MPA
AllocateBoats <- function(nboats, IsFished, pop)  {
   FishableBiomass <- pop*IsFished     # get relative fishable biomass
   ## how many boats to allocate to each area per unit FishableBiomass
   scale <- nboats/sum(FishableBiomass)
   boats <- FishableBiomass * scale
   return(boats)
}
temp <- IsFished(ncells=8, MPA.width=3)
x <- AllocateBoats(nboats=100, IsFished=temp, pop=rep(100,8))



#Add boats to MPA model with logistic growth
#r=rate of increase, K=carrying capacity, 
#ncells=number of cells for model, nsteps=number of time periods,
#MPA.width=number of cells in MPA, mrate=movement rate
MPA.model.lab8.boats <- function(r, K, u.out, ncells, nyear, MPA.width, sigma, nboats) {
   #fraction of fish that each boat will catch
   boatQ <- u.out*ncells/nboats  
   is.fished <- IsFished(ncells=ncells, MPA.width=MPA.width)
   
   #create a matrix N with nyear rows and ncells columns
   N <- matrix(data=0, nrow=nyear, ncol=ncells)
   N[] <- 0
   
   #set the abundance in the first year to K in each cell
   N[1,] <- K
   
   #calculate the M matrix here (only need do this once!)
   M <- DispersalKernel(sigma=sigma, ncells=ncells)
   
   #plot the initial numbers
   plot(x=1:ncells, y=N[1,], xlab="Cell number", lwd=3,
        ylab="Population size", ylim=c(0, 1.05*max(N)), type="l", yaxs="i", xaxs="i")
   
   for (i in 1:(nyear-1)) {
      #calculate surplus production from the logistic model
      surplus.prod <- r*N[i,]*(1-N[i,]/K)
      
      #calculate catch numbers from harvest rate
      boats <- AllocateBoats(nboats=nboats, IsFished=is.fished, pop=N[i,])
      u.vec <- 1-exp(-boats*boatQ)
      catches <- u.vec*N[i,]
      
      #calculate pre-moving abundance in each cell (store in temp vector)
      preN <- N[i,] + surplus.prod - catches
      
      #use matrix multiplication to disperse premove abundance to cells next year
      N[i+1,] <- preN %*% M
      
      lines(x=1:ncells, y=N[i+1,], lwd=(nyear-i+1)/nyear*3)
   }
}

MPA.model.lab8.boats(r=0.2, K=1000, u.out=0.3, ncells=21, MPA.width=5, 
                     nyear=10, sigma=0.5, nboats=100)
