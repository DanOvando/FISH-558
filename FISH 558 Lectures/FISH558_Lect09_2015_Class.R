###############################################################################
#### FISH 558
#### Lecture 09
#### Kelli Faye Johnson
###############################################################################

# =============================================================================
# Step 01: add the functions to your environment
# =============================================================================

# Slide 04
# Population model with basic specifications
#' @param r A numeric value specifying the population growth rate, where the
#'   population experiences exponential growth.
#' @param Start_val A numeric value specifying the initial population size.
#' @param Thresh A numeric value specifying the quasi-extinction level, or
#'   population size at which quasi-extinction occurs.
#' @param Nyears A numeric value specifying the number of years in a single
#'   simulation.
#' @param Nsim A numeric value specifying the number of simulations to perform.
#' @param SigP A numeric value specifying the yearly variation in population.
DoProject <- function(r,Start_val=500,Thresh=10,Nyears=100,Nsim=1000,SigP=0.2)
{
 WentExtinct <- 0
 for (Isim in 1:Nsim)
  {
   Nvalue <- Start_val
   Extinct <- F
   for (Iyr in 1:Nyears)
    {
     Nvalue <- r*Nvalue*exp(rnorm(1,0,SigP) - SigP^2/2)
     if (Nvalue < Thresh) Extinct <- T
    }
   if (Extinct==T) WentExtinct <- WentExtinct + 1
  }
 return(WentExtinct/Nsim*100)

}

# Slide 17
# Iberian lynx function
#' @param K
#' @param s24 A vector of numeric values specifying the survival rate of ages two
#'   through four, thus the vector has a length of three. There are no default
#'   values.
#' @param s1 A vector of two numeric values specifying the survival rate of cubs.
#'   The first entry is the survival rate of cubs in the absence of a catastrophe
#'   and the second value is the suvival rate of cubs during a catastrophe.
#'   There are no default values.
#' @param Ninit
#' @param ProbFemCub A numeric value specifying the probability that a cub is female.
#'   The default value is \code{0.87}.
#' @param y A numeric value between zero and one, specifying the probability of
#'   a catastrophe, where catastrophic events are caused by droughts.
#' @param Thresh An integer specifying the threshold level, below which quasi-extinction
#'   occurs. The default value is \code{2}.
#' @param Nsim A numeric value specifying the number of simulations.
#' @param Nyear A numeric value specifying the number of years in a simulation.
Lynx <- function(K, s24, s1, Ninit, ProbFemCub = 0.87, ProbCat = 0.1,
  Thresh = 2, Nsim = 100,Nyear = 50) {
 WentExtinct <- 0
 for (Isim in 1:Nsim) {
    # Initialize
    n <- Ninit
    Extinct <- F
    for (Iyr in 1:Nyear) {
      # Define year-1 survival
      if (runif(1,0,1) < ProbCat) SurvJ <- s1[2] else SurvJ <- s1[1]
      s <- c(SurvJ,s24)

      # Update the stage-structure
      Nnext <- rep(0,4)
      for (Istage in 1:4)
       Nnext[Istage] <- rbinom(1,n[Istage],s[Istage])

      # Generate the cub
      Cubs <- if (n[4] > 0) rbinom(1,n[4],ProbFemCub) else 0

      # Probaby becoming a breeder
      if (Nnext[2]+Nnext[3] > 0) ProbBreed <- (K-Nnext[4])/(Nnext[2]+Nnext[3]) else ProbBreed <- 0.5
      if (ProbBreed < 0) ProbBreed <- 0
      if (ProbBreed > 1) ProbBreed <- 1
      Ntemp <- Nnext[2]+Nnext[3]
      NewBreed <- if (Ntemp > 0) rbinom(1,Ntemp,ProbBreed) else 0
      n[4] <- Nnext[4] + NewBreed
      n[3] <- Nnext[2] + Nnext[3] - NewBreed
      n[2] <- Nnext[1]
      n[1] <- Cubs

      # Check for extinction
      Ntot <- sum(n)
      if (Ntot <= Thresh) Extinct <- T
    }
    if (Extinct==T) WentExtinct <- WentExtinct + 1
 }
 return(WentExtinct/Nsim*100)

}

# =============================================================================
# Step 02: Run the initial model
# =============================================================================
# x is a vector of intrinsic growth rates
x <- c(0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15)
# Initiate the probability of extinction to zero for each growth rate
y <- rep(0, length = length(x))
for (II in 1:length(x)) {
  y[II] <- DoProject(x[II])
}
 cat("The probability of extinction given the specified growth rate.\n")
 print(y)

 plot(x, y, xlab = "r", ylab = "Probability (extinction)", las = 1,
   type = "b", lty = 1, lwd = 2, pch = 16)

# =============================================================================
# Step 03: Change DoProject to allow for temporal correlation
# =============================================================================

 # Slide 04
 # Population model with basic specifications
 #' @param r A numeric value specifying the population growth rate, where the
 #'   population experiences exponential growth.
 #' @param Start_val A numeric value specifying the initial population size.
 #' @param Thresh A numeric value specifying the quasi-extinction level, or
 #'   population size at which quasi-extinction occurs.
 #' @param Nyears A numeric value specifying the number of years in a single
 #'   simulation.
 #' @param Nsim A numeric value specifying the number of simulations to perform.
 #' @param SigP A numeric value specifying the yearly variation in population.
 DoAutoProject <- function(r,Start_val=500,Thresh=10,Nyears=100,Nsim=1000,SigP=0.2,p = .717,catastrophe = 0)
 {
   WentExtinct <- 0

   for (Isim in 1:Nsim)
   {
     Nvalue <- Start_val

      Extinct <- F

     errors <- matrix(0,1,Nyears +1)

     for (Iyr in 1:Nyears)
     {
       errors[Iyr+1] <- p*errors[Iyr] + sqrt(1-p^2)*rnorm(1,0,SigP)
       Nvalue <- r*Nvalue*exp(errors[Iyr+1] - SigP^2/2)

        if (Nvalue < Thresh) Extinct <- T
     }
     if (Extinct==T) WentExtinct <- WentExtinct + 1
   }
   return(WentExtinct/Nsim*100)

 }


 #With autocorrelation

 # x is a vector of intrinsic growth rates
 x2 <- c(0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15)
 # Initiate the probability of extinction to zero for each growth rate
 y2 <- rep(0, length = length(x))
 for (II in 1:length(x)) {
   y2[II] <- DoAutoProject(x2[II])
 }
 cat("The probability of extinction given the specified growth rate.\n")
 print(y2)


 with.ar <- data.frame(scenario = 'with AR', x = x2, y = y2)

 without.ar <- data.frame(scenario = 'without AR', x = x, y = y)

 results <- rbind(with.ar,without.ar)

 ggplot(results,aes(x,y/100,color = scenario, shape = scenario)) +
   geom_line(aes(linetype = scenario),size = 2, alpha = 0.6) +
   geom_point(size = 4) +
   xlab('r') +
   ylab("Probability (extinction)") +
   scale_y_continuous(labels = percent)


quartz()
plot(x2, y2, xlab = "r", ylab = "Probability (extinction)", las = 1,
      type = "b", lty = 1, lwd = 2, pch = 16)
par(new = T)
plot(x, y, xlab = "r", ylab = "Probability (extinction)", las = 1,
      type = "b", lty = 2, lwd = 2, pch = 16)

# =============================================================================
# Step 04: Run the Iberian lynx example
# =============================================================================
# xx is a vector of hypothesized number of territories
xx <- 1:10
# Initiate the probability of extinction at zero for each number of territories
yy <- rep(0, length = length(xx))
for (II in 1:length(xx)) {
  yy[II] <- Lynx(
    K = xx[II],
    s24 = c(0.7,0.7,0.86),
    s1 = c(0.5,0.2),
    Ninit = c(4,2,0,5),
    ProbFemCub = 0.87,
    ProbCat = 0.1,
    Thresh = 2,
    Nsim = 1000,
    Nyear = 50
  )
}

cat("The probability of extinction given the specified number of territories.\n")
print(yy)

 plot(xx, yy, xlab = "Number of territories", lty = 1, lwd = 2, pch = 16,
 type = "b", ylab = "Probability (extinction)", las = 1)
