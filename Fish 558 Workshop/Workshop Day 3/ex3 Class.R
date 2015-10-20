# Ex3 <- function()
#

#   # call for the logistic model
  set.seed(443)
  SchModel <- DoSir(Nout=100,Model='Schaefer')

  head(SchModel)

  quartz()
  ggplot(SchModel,aes(K,r)) + geom_hex()

  ggplot(SchModel,aes(x = K,y = r)) + geom_point(aes(fill = NegLogLike), shape = 21, size = 2) + scale_fill_gradient(low = 'green',high = 'red')

  ggplot(SchModel,aes(NegLogLike)) + geom_histogram()

  # call for the exponential model
  ExpoModel <- DoSir(Nout=100,Model= 'ExponModel')



# =================================================================================

DoSir <- function(Nout=1000,Model)
{

  # Read in the basic data
  TheData <- ReadData()
  Yr1 <- TheData$CatchYr[1]
  Catch <- TheData$CatchVal
  SurveyEst <- TheData$SurveyEst
  SurveyCV <- TheData$SurveyCV
  SurveyYr <- TheData$SurveyYr
  Nyears <- length(Catch)
  years <- TheData$CatchYr

  # Storage for the parameters and the final depletion
  Vals <- as.data.frame(matrix(0,ncol=5,nrow=Nout))

  colnames(Vals) <- c('K','r','Pop1965','AddCV','NegLogLike')
  # Storage for the total likelihood encountered in the first stage sampling
  # and the number of draws
  AveLike <- 0
  Ntest <- 0

  # Reset parameters for SIR
  Threshold <- exp(0)
  Cumu <- 0
  Ndone <- 0
  while (Ndone < Nout)
  {
    # Generate from priors

    r <- runif(1,0,.15)

    Pop1965 <- runif(1,10000,15000)

    AddCV <- runif(1,.1,.2)

    K <- runif(1,20000,50000)

    # Call the population model
    Pop <- PopModel(Catch = Catch,r = r,K = K,years = years,InitPop = Pop1965,ExponModel = Model)

    survey <- data.frame(TheData$SurveyYr,TheData$SurveyEst,TheData$SurveyCV)

    colnames(survey) <- c('year','SurveyEst','SurveyCV')

    Pop <- join(Pop,survey, by = 'year')

    #     ggplot(Pop,aes(year,n)) + geom_point() + geom_line(aes(year,SurveyEst))
    #

    #write the pop model
    # Compute the negative log-likelihood and hence the likelihood
    NegLogLike <- Likelihood(Pop = Pop$n,SurveyYr-Yr1+1,SurveyEst,SurveyCV,AddCV)
    TheLike <- exp(-1*NegLogLike-32.19)

    # Determine if a parameter vector is to be saved
    Cumu <- Cumu + TheLike

    AveLike <- AveLike + TheLike
    Ntest <- Ntest +1

    while (Cumu > Threshold & Ndone < Nout)
    {
      Ndone <- Ndone + 1
      Cumu <- Cumu - Threshold

      Vals[Ndone,] <- data.frame(K,r,Pop1965,AddCV,NegLogLike)

    }
  }

  Vals$AveLike <- AveLike/Ntest

  return(Vals)
}

# =================================================================================

Likelihood <- function(Pop,SurveyYr,SurveyEst,SurveyCV,AddCV)
{
  # Account for the additional CV
  UseCV <- sqrt(SurveyCV^2+AddCV^2)

  # Extract the predictions corresponding to the observations and compute the negatuve log-likelihood
  Preds <- Pop[SurveyYr]
  Residuals <- log(UseCV)+0.5*(log(Preds)-log(SurveyEst))^2/UseCV^2
  LogLike <- sum(Residuals)
}

# =================================================================================
PopModel <- function(Catch,r,K,years,InitPop,ExponModel)
{

  time <- length(years)

  output <- as.data.frame(matrix(NA,nrow = (time),ncol = 3))

  colnames(output) <- c('year','catch','n')

  output$catch <- Catch

  output$n[1] <- InitPop

  output$year <- years

  if (ExponModel == 'ExponModel')
  {

    for (t in 2:time)
    {

      output$n[t] <- pmax(1e-5,(1+r)*output$n[t-1] - output$catch[t-1])

    }

  }
  if (ExponModel == 'Schaefer')
  {

    for (t in 2:time)
    {
      output$n[t] <- pmax(1e-5,output$n[t-1] + (output$n[t-1]*r)*(1-output$n[t-1]/K) - output$catch[t-1])
    }

  }

  return(output)
}

# =================================================================================
ReadData <- function()
{
  TheData1 <- read.csv(paste('Fish 558 Workshop/',lecture,'/Ex3a.csv', sep = ''),header=TRUE, stringsAsFactors = F)

  TheData2 <- read.csv(paste('Fish 558 Workshop/',lecture,'/Ex3b.csv', sep = ''),header=TRUE, stringsAsFactors = F)

  Outs <- NULL
  Outs$SurveyYr <- TheData1[,1]
  Outs$SurveyEst <- TheData1[,2]
  Outs$SurveyCV <- TheData1[,3]
  Outs$CatchYr <- TheData2[,1]
  Outs$CatchVal <- TheData2[,2]
  return(Outs)
}

# =================================================================================
Ex3()
