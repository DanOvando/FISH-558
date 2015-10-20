DoSir <- function(Nout=1000,dat,r.manual = NA, K.manual = NA, phi.manual = NA,progbar = F)
{
  if (progbar == T)
  {
    pb <- txtProgressBar(min = 1, max = Nout, style = 3)
  }
  # Read in the basic data
  Yr1 <- dat$year[1]
  Catch <- dat$catch
  SurveyEst <- dat$catch.rate
  #   SurveyCV <- TheData$SurveyCV
  #   SurveyYr <- TheData$SurveyYr
  Nyears <- length(Catch)
  years <- dat$year

  # Storage for the parameters and the final depletion
  Vals <- as.data.frame(matrix(0,ncol=6,nrow=Nout))

  colnames(Vals) <- c('K','r','phi','sigma','biomass','thelike')
  # Storage for the total likelihood encountered in the first stage sampling
  # and the number of draws
  AveLike <- 0
  Ntest <- 0

  # Reset parameters for SIR
  Threshold <- exp(-77)
  Cumu <- 0
  Ndone <- 0
  while (Ndone < Nout)
  {
    # Generate from priors
    if (is.na(r.manual))
    {
      r <- runif(1,.15,.45)
    }
    if (is.na(phi.manual))
    {
      phi <- runif(1,0.8,1.2)
    }
    if (class(phi.manual) == 'numeric')
    {
      phi <- phi.manual
    }
    #     AddCV <- runif(1,.1,.2)
    if (is.na(K.manual))
    {
      K <- runif(1,5000,8000)
    }
    # Call the population model
    #     browser()

    #     optim(par = c(log(.3),log(5000),log(1)),schaefer.SIR.model,dat = dat, time = time)
    Pop <- schaefer.model(pars = log(c(r,K,phi)),dat = dat,time = years, mode = 'popmodel')

    TheLike <- Pop$total_like[1]  #exp(-1*NegLogLike-32.19)

    Cumu <- Cumu + TheLike

    AveLike <- AveLike + TheLike

    Ntest <- Ntest +1
    while (Cumu > Threshold & Ndone < Nout) #check and see if cumulative likelihood is over threshold
    {
      Ndone <- Ndone + 1

      if (progbar == T)
      {
      setTxtProgressBar(pb, Ndone)
      }

      Cumu <- Cumu - Threshold
      Vals[Ndone,] <- data.frame(K,r,phi,Pop$sigma,last(Pop$biomass),TheLike)

    }
  }

  Vals$AveLike <- AveLike/Ntest

  Vals$BvK <- Vals$biomass/Vals$K
  return(Vals)
}