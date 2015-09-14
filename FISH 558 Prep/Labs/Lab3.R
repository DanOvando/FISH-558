## ---------- FISH 458 | SPRING 2012 | University of Washington ----------
## Lab 3 : Likelihood profiles
## Cole Monnahan | monnahc@uw.edu
## 4/10/2012

## This lab examines CPUE of a Hake fishery and shows how to
## perform likelihood profile calculations in R.

## The source() function loads in functions saved in an external .R
## file. This particular file contains an edited version of the
## filled.contour function that removes the Z-scale, which screws up
## plotting points on top of the surface. This file needs to be in the
## workspace.
source("my.filled.contour.R")
library(stats4)

## Load in the data and create vectors from it for easier use. I've
## capitalized them to distinguish externally load data from something
## I create in R, or function parameters.
Data <- read.csv("hake_data.csv")
Years <- Data[,1]
Catches <- Data[,2]
CPUE.obs <- Data[,3]

## ------------------------------------------------------------

## ---------- Begin functions that need to be run once
## First, the  logistic growth model that calculates CPUE.

CPUE.model <- function(r, K, q, year.start=1965, year.end=1987,
                       catches){
    ## This function returns a vector of biomass (B) and CPUE (q*B) for
    ## the years given

    ## First setup the variables needed
    years <- year.start:year.end
    num.years <- length(years)
    biomass <- rep(NA, num.years)      # vector of abundances
    cpue <- rep(NA, num.years)      # vector of CPUE

    ## Now calculate the population and CPUE
    biomass[1] <- K
    cpue[1] <- q*K
    for(i in 2:num.years){
        Bt <- biomass[i-1]
        ## Use max() to keep it from going negative, .01 is arbitrary
        biomass[i] <- max(.01,Bt+Bt*r*(1-Bt/K)-catches[i])
        cpue[i] <- biomass[i]*q
    }

    ## Return the results as a matrix of 2 columns
    result <- cbind(biomass, cpue)

    return(result)
}

CPUE.mle <- function(r.log, K.log, q.log,  catches, CPUE){
    ## This function is only used to find the optimal parameter
    ## values. Note that the parameters for this function are in
    ## log-space, which is convenient for optimization when parameters
    ## should be positive. It also puts them on a much closer scale
    ## which helps convergence.

    ## Undo the log since our logistic model needs them to be.
    r <- exp(r.log)
    K <- exp(K.log)
    q <- exp(q.log)

    ## Run the model and grab predicted CPUE
    CPUE.temp <- CPUE.model(r=r, K=K, q=q, catches=catches)[,2]

    ## Calculate SSQ
    SSQ <- sum( ( log(CPUE.obs)-log(CPUE.temp))^2 )

    return(SSQ)
}
## ---------- End functions



## Now we can explore the models. We can create a "surface" in 2D by
## calculating SSQ at a whole range of (r,K, q)
## coordinates. Unfortunately this means we can only vary two
## parameters. From playing around (in Excel or in R) I know that
## .0005 is a reasonable value for q so I'll fix it at that and see
## how r and K affect SSQ.

## Setup a grid, just like in Excel.
num.seq <- 100
r.seq1 <- seq(0,1.5, len=num.seq)
K.seq1 <- seq(0, 7000, len=num.seq)
SSQ.mat <- matrix(NA, nrow=num.seq, ncol=num.seq)

## Now a double loop runs through each cell in the matrix and
## calculates SSQ and stores it.
for(i in 1:num.seq){
    r.temp <- r.seq1[i]
    for(j in 1:num.seq){
        K.temp <- K.seq1[j]
        SSQ.mat[i,j] <- CPUE.mle(r.log=log(r.temp), K.log=log(K.temp),
                                 q.log=log(.0005), catches=Catches,
                                 CPUE=CPUE.obs)
    }
}

## Now scale it and plot surface, note that I take the log of SSQ to
## highlight the differences in SSQ.
SSQ.max <- max(SSQ.mat, na.rm=TRUE)
SSQ.mat <- SSQ.mat/SSQ.max
my.filled.contour(x=r.seq1, y=K.seq1, z=log(SSQ.mat),
                  main="Likelihood Surface for r,K when q=.0005",
                  xlab="r", ylab="K")
## In this plot the blues are "good" fits and the pinks are really
## bad. You can see that there are regions of blue that are pretty
## close to the same color. This means that the SSQ would be about the
## same in this region, meaning a range of models fit about as well.

## Now let's use mle() to find the minimum of this surface. These
## following vectors are my "test" starting points. They should all
## converge to the same point (the global minimum).
r.test <- c(.1,.2, .4, .8, 1,1.5, .6, .3,1)
K.test <- c(7,5,4,3,2,1,6,6,5)*1000
q.test <-  c(4, 3,5, 6, 10, 1, 3, 5, 4)/10000
test.mat <- matrix(NA, nrow=length(r.test), ncol=4)

## Run the optimizer for each of these parameter combinations
for(i in 1:length(r.test)){
    MLE.test <- mle(CPUE.mle, start=list(r.log=log(r.test[i]),
                         K.log=log(K.test[i]), q.log=log(q.test[i])),
               fixed=list(catches=Catches,
               CPUE=CPUE.obs))

    test.mat[i,1] <- exp(coef(MLE.test)[1])
    test.mat[i,2] <- exp(coef(MLE.test)[2])
    test.mat[i,3] <- exp(coef(MLE.test)[3])
    test.mat[i,4] <- -as.numeric(logLik(MLE.test))
}

test.mat
## Yes they all converge except one, which is obviously not as
## good. We can tell because the SSQ (the last column) is much higher.


## Visual proof:
points(r.test, K.test, pch=16)          # Starting values
arrows(x0=r.test, y0=K.test, x1=test.mat[,1], y1=test.mat[,2],
       length=.1)                       # arrows show movement



## Let's start it again from the MLE from above just to make sure it
## doesn't move, and then use this as a final MLE

MLE <- mle(CPUE.mle, start=list(r.log=log(test.mat[1,1]),
                          K.log=log(test.mat[1,2]), q.log=log(test.mat[1,3])),
                fixed=list(catches=Catches,
                CPUE=CPUE.obs))

## The optimal parameters
(r.hat <- as.numeric(exp(coef(MLE)[1])))
(K.hat <- as.numeric(exp(coef(MLE)[2])))
(q.hat <- as.numeric(exp(coef(MLE)[3])))
(SSQ.min <- -as.numeric(logLik(MLE)))   # The minimum SSQ

## Let's see what it looks like compared to the data
plot(Years, CPUE.obs, pch=16, ylim=c(0,2),
     main="Hake CPUE: Observed vs. MLE", ylab="CPUE")
lines(Years,CPUE.model(r=r.hat, K=K.hat, q=q.hat, catches=Catches)
      [,2], lwd=2, col="red")
legend("topright", legend=c("Observed CPUE", "Modeled CPUE"),
       lty=c(NA, 1), col=c(1,2), pch=c(16, NA), lwd=2)

## ---------- PROFILE of r
num.seq <- 50
r.seq <- seq(.01, 1.5, len=num.seq)
r.profile <- matrix(NA, nrow=num.seq, ncol=4)
r.profile[,1] <- r.seq

## Lets run it once to make sure our first value converges, and so
## that we can use these as our starting seeds for the second step
mle.temp <-mle(CPUE.mle, start=list( K.log=log(7500),
                q.log=log(q.hat)), fixed=list(r.log=log(r.seq[1]),
                catches=Catches, CPUE=CPUE.obs))
(r.profile[1,1] <- as.numeric(exp(coef(mle.temp)[1])))
(r.profile[1,2] <- as.numeric(exp(coef(mle.temp)[2])))
(r.profile[1,3] <- as.numeric(exp(coef(mle.temp)[3])))
(r.profile[1,4] <- as.numeric(-logLik(mle.temp)))

## Now loop through each value of r and recalculate the other
## parameters that minimize SSQ
for(i in 2:num.seq){
    r.temp <- r.seq[i]
    mle.temp <- NULL                    # reset this
    ## Start the next value of r from the previous values of K and q,
    ## which should be pretty close.
    mle.temp <- mle(CPUE.mle, start=list(K.log=log(r.profile[i-1,2]),
                              q.log=log(r.profile[i-1,3])),
                              fixed=list(r.log=log(r.temp),
                              catches=Catches, CPUE=CPUE.obs))
    r.profile[i,4]<- -as.numeric(logLik(mle.temp))

    r.profile[i,2]<- exp(coef(mle.temp)[2])
    r.profile[i,3] <- exp(coef(mle.temp)[3])

}

## Plot the results
par(mfrow=c(2,2),oma=c(0,0,1.5,0))
plot(r.seq, r.profile[,4], type="l", lwd=2,xlab="r",ylab="SSQ",
     main="r vs. SSQ")
abline(v=r.hat, col="red")              # This was the min from before
plot(r.seq, r.profile[,2], type="l", lwd=2, xlab="r",
     ylab="K",main="r vs. K")
plot(r.seq, r.profile[,3], type="l", lwd=2, xlab="r",
     ylab="q",main="r vs. q")
plot(r.profile[,2], r.profile[,3], type="l", lwd=2, xlab="K",
     ylab="q",main="K vs. q")
mtext(text="r Profile", outer=TRUE, line=-1, cex=1.5)

## What does this look like on the likelihood surface?
par(mfrow=c(1,1))
my.filled.contour(x=r.seq1, y=K.seq1, z=log(SSQ.mat),
                  main="Likelihood Surface for r,K when q=.0005",
                  xlab="r", ylab="K")
points(r.seq, r.profile[,2])
## Why doesn't this line follow the contours exactly? Should it? Think
## about how I created this surface.

## The other parameter profiles can be found in a similar fashion.

## End of lab
## ------------------------------------------------------------
