#####################################################################
#Homework 4: estimating MSY, Bmsy, umsy from a 
#general stock assessment model
#FISH 458 School of Aquatic and Fishery Sciences
#University of Washington
#Your name, your email address
#Dates your worked on this
#####################################################################

#list of model parameter values
alpha <- 0.0001   #length-weight alpha
beta <- 3.05      #length-weight beta
Linfinity <- 80   #von Bertalanffy Linfinity
k <- 0.2          #von Bertalanffy K
t0 <- -0.2        #von Bertalanffy t0
h <- 0.7          #Beverton-Holt steepness is R/R0 when S = 0.2SSB0. 
R0 <- 1000000     #number of recruits in unfished population
Lv50 <- 40        #length at which 50% of fish are vulnerable to gear
Lv95 <- 45        #length at which 95% of fish are vulnerable to gear
pfemale <- 0.5    #proportion of the population that is female
Lm50 <- 30        #length at which 50% of females are mature
Lm95 <- 35        #length at which 95% of females are mature
surv <- 0.8       #natural survival, assumed constant for all ages
n <- 10           #plus group age 

#=========PART ONE================================================================
#Calculates unfished equilibrium numbers at age, including the plus group.
#Inputs: unfished rec, natural survival, and
#plus group age. Assumes R0 is for age 1 individuals. 
#=================================================================================
unfished.num <- function(R0, surv, n) {
   N.vec <- vector(length=n)

   
   
   return(N.vec)
}
x <- unfished.num(R0=R0, surv=surv, n=n)
print(x)

#=========PART TWO================================================================
#given a vector of ages, and parameters, returns a vector of weights at those 
#ages. 
#Inputs: age.vec a vector of ages
# alpha, beta of length-weight; Linfinity, k, t0 of Von Bertalanffy
#=================================================================================
weight.at.age <- function(age.vec, alpha, beta, Linfinity, k, t0) {
   N.ages <- length(age.vec)

   
   
   return(weight.vec)
}
x <- weight.at.age(age.vec=1:n,alpha=alpha, beta=beta, Linfinity=Linfinity, k=k,t0=t0)
print(x)

#=========PART THREE==============================================================
#Proportion of population at age that is spawning females
#Given a vector of ages, returns the proportion at each that is female and spawning.
#Assumes a logistic curve with Lm50 the length at which 50% are mature and 
#Lm95 the length at which 95% are mature
#Inputs: age.vec a vector of ages, Linfinity, k, t0, Lm50, Lm95, pfemale
#=================================================================================
prop.fem.spawners <- function(age.vec, Linfinity, k, t0, Lm50, Lm95, pfemale) {

   
   
   return(prop.mature)
}
x <- prop.fem.spawners(age.vec=1:n,Linfinity=Linfinity, k=k,t0=t0, Lm50=Lm50, 
                       Lm95=Lm95, pfemale=pfemale)
print(x)

#=========PART FOUR===============================================================
#Calculate SSB0, which comes from the sum of the spawning biomass per recruit
#multiplied by the number of recruits. In this model this is also the unfished 
#numbers-at-age multiplied by the weight-at-age 
#multiplied by the proportion of the fish that are spawning females.
#Inputs: n, alpha, beta, Linfinity, k, t0, Lm50, Lm95, pfemale, surv, R0
#=================================================================================
calc.SSB0 <- function(n, alpha, beta, Linfinity, k, t0, Lm50, Lm95, pfemale, surv, R0) {

   
   
   return(SSB0)
}
x <- calc.SSB0(n=n, alpha=alpha, beta=beta, Linfinity=Linfinity, k=k, t0=t0, 
          Lm50=Lm50, Lm95=Lm95, pfemale=pfemale, surv=surv, R0=R0)
print(x)

#==========PART FIVE===============================================================
#Vulnerability to fishing, as a function of length, converts an age-vector into 
#vulnerability values by age. First converts age to length, then converts length
#to vulnerability assuming a logistic selectivity curve 
#=================================================================================
prop.vulnerable <- function(age.vec, Linfinity, k, t0, Lv50, Lv95) {

   
   
   return(prop.vuln)
}
x <- prop.vulnerable(age.vec=1:n,Linfinity=Linfinity, k=k,t0=t0, Lv50=Lv50, Lv95=Lv95)
print(x)

#==========PART SIX===============================================================
#Numbers of fish in the population as function of harvest rate u, recruits R, plus 
#age a, the vulnerability parameters, and von Bertalanffy parameters.
#Returns a vector of numbers at age
#===================================================================================
calc.pop.size <- function(R, n, u, surv, Linfinity, k, t0, Lv50, Lv95) {

   
   
   return(N.vec)
}
x <- calc.pop.size(R=1, n=n, u=0, surv=surv, Linfinity=Linfinity, k=k, 
              t0=t0, Lv50=Lv50, Lv95=Lv95)
print(x)

#==========PART SEVEN==============================================================
#Calculate equilibrium yield and spawning biomass. 
#Requires all the parameters of the model
#Calls all the little bits and pieces calculated so far. 
#==================================================================================
eqm.Y.SSB <- function(alpha, beta, Linfinity, k, t0, h, R0, Lv50, Lv95,
                     pfemale, Lm50, Lm95, surv, n, u) {
   
   #numbers at age from 1 recruit as a function of harvest rate

   
   #proportion of spawners, and weight at age

   
   #SBPR(u) spawning biomass per recruit as a function of u
   
   
   #YPR(u) yield per recruit at harvest rate u
   
   
   #caculate SSB0 (for Beverton-Holt)
   
   
   #calculate equilibrium number of recruits from Beverton-Holt and SBPR(u)
   
   
   #calculate total equilibrium yield and spawning biomass at this exploitation rate
   
   
   #returns a list, to access elements use x[[1]] and x[[2]]
   return(list(yield=yield, SSB=SSB, YPR=YPR, SBPR=SBPR, eqm.recruits=eqm.recruits))
   
}
#calculate equilibrium yield and spawning biomass for a given harvest rate
x <- eqm.Y.SSB(alpha=alpha, beta=beta, Linfinity=Linfinity, k=k, t0=t0,
                    h=h, R0=R0, Lv50=Lv50, Lv95=Lv95,
                    pfemale=pfemale, Lm50=Lm50, Lm95=Lm95, surv=surv, n=n, u=0.1)
print(x)

#==========PART EIGHT==============================================================
#Calculate yield, SSB for a range of exploitation rates from 0 to 1
#Requires all the parameters of the model
#Calls all the little bits and pieces calculated so far. 
#==================================================================================
calc.ref.points <- function(alpha, beta, Linfinity, k, t0, h, R0, Lv50, Lv95,
                             pfemale, Lm50, Lm95, surv, n) {
   
   
   
   results <- matrix(nrow=ncalcs, ncol=3)  #ncalcs is number of rows
      
   
   return(results)
}
x <- calc.ref.points(alpha=alpha, beta=beta, Linfinity=Linfinity, k=k, t0=t0,
                           h=h, R0=R0, Lv50=Lv50, Lv95=Lv95,
                           pfemale=pfemale, Lm50=Lm50, Lm95=Lm95, surv=surv, n=n)
head(x)   #print out the first few rows
#write the results to a file for later inspection to find the MSY etc. 
write.csv(x=x, file="MSY results.csv")   
   