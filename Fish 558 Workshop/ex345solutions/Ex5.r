Project <- function(Props,Level,States,NumInEachState,Print)
{  
 Nstate <- length(NumInEachState)
 Nprop <- length(Props) 
  
 ExpectedValue <- 0
 TableLine <- rep(0,Nstate)
 
 # Proportion which rrecovert
 Recover = 1.0/(1+exp(-1*0.2*(Level-10)))
 
 # Work through each entry
 for (II in 1:Nprop)
  {
   # Number of people without the disease
   Perf <- (1-Props[II]) + Props[II] * Recover
   # Objective function
   Perf <- Perf*10000 - Level*100
   TableLine[States[II]] <- TableLine[States[II]] + Perf
   ExpectedValue <- ExpectedValue + Perf
 } 
 
 # Find expected values
 ExpectedValue <- ExpectedValue / Nprop  
 for (JJ in 1:Nstate)
  TableLine[JJ]  <- TableLine[JJ] / NumInEachState[JJ]
 if (Print==T) cat(Level,TableLine,ExpectedValue,"\n")
 return(ExpectedValue)
 
}

xx <- rbeta(1000000,12,36)

yy <- dbinom(45,150,xx)

zz <- sample(xx,500000,replace=T,prob=yy)

par(mfrow=c(2,2))
hist(xx)
hist(zz)

# Part A - check the posterior
alpha <- 12+45
beta <- 36+(150-45)
meanV <- alpha/(alpha+beta)
var <- alpha*beta/((alpha+beta)^2*(alpha+beta+1))
cat(alpha,beta,meanV,sqrt(var),"\n")
cat("Check",mean(zz),sd(zz),"\n")

# Select the cut-offs
pers <- c(0.2,0.4,0.6,0.8)
Cutoff <- qbeta(pers,alpha,beta)
Ncut <- length(Cutoff)
print(Cutoff)

# Generate from the posterior
Props <- rbeta(2000,alpha,beta)
hist(Props)
Nprop <- length(Props)

# Assign each proportion to a state of nature
States <- rep(0,Nprop)
NumInEachState <- rep(0,Ncut+1)
for (II in 1:Nprop)
 { 
  Icol <- Ncut + 1
  for (JJ in Ncut:1) if (Props[II] < Cutoff[JJ]) Icol <- JJ  
  States[II] <- Icol
  NumInEachState[Icol] <- NumInEachState[Icol] + 1
 }
print(NumInEachState)

Project(Props,4,States,NumInEachState,Print=T)
Project(Props,10,States,NumInEachState,Print=T)
Project(Props,14,States,NumInEachState,Print=T)
Project(Props,20,States,NumInEachState,Print=T)

Vals <- seq(from=0,to=20,by=1)
Perf <- rep(0,length(Vals))
for (II in 1:length(Vals))
 Perf[II] <-  Project(Props,Vals[II],States,NumInEachState,Print=F) 
plot(Vals,Perf,xlab="Level of Mitigation",ylab="Performance metric",type='b')


