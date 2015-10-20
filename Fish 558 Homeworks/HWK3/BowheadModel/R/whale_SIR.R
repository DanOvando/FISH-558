#' Do a simple SIR to test validity of parameters
#'
#' \code{do.simple.sir} returns a posterior of calf and adult
#' survival, max fecundity, and K, such that the population
#' does not crash given an observed catch history
#' @param Nout the number of desired draws from the posterior
#' @param dat the raw whale data
#' @param s.0.lower lower bound of calf survival
#' @param s.1.lower upper bound of calf survival
#' @param mode 'simple' or 'fit' tells it to either run
#' SIR on crashing th population, 'fit' tells it to fit the
#' model
#' @export

whale.sir <- function(Nout=1000,dat,catch.dat = catch.dat,s.0.lower = 0.8,s.0.upper = 1,
                      s.rest.lower = 0.9, s.rest.upper = 1,
                      k.lower = 10000,k.upper = 20000,
                      f.max.lower = 0.25,f.max.upper = 0.33,progbar = F,mode = 'simple', MLE = 2.22)
{
  if (progbar == T)
  {
    pb <- txtProgressBar(min = 1, max = Nout, style = 3)
  }


  # Storage for the parameters and the final depletion
  Vals <- as.data.frame(matrix(0,ncol=6,nrow=Nout))

  colnames(Vals) <- c('s.0','s.rest','K','f.max','whales','thelike')
  # Storage for the total likelihood encountered in the first stage sampling
  # and the number of draws
  AveLike <- 0
  Ntest <- 0

  Threshold <- exp(0)
  Cumu <- 0
  Ndone <- 0
  while (Ndone < Nout)
  {
    # Generate from priors

    s.0.guess <- runif(1,s.0.lower, s.0.upper)

    s.rest.guess <- runif(1,s.rest.lower, s.rest.upper)

    K.guess <- runif(1,k.lower, k.upper)

    f.max.guess <- runif(1,f.max.lower,f.max.upper)

    whales <- make.whales(dat = dat,catch.dat = catch.dat,f.max = f.max.guess,s.0 = s.0.guess,
                          s.rest = s.rest.guess, K = K.guess)

    #     whales <- make.whales(dat = dat,catch.dat = catch.dat,f.max = .26,s.0 = .84,s.rest = .8, extra.time = 0,Kinit = 1846)

    #Run the model using the created object
    whale.proj <- whale.pop.model(whales = whales, dat = dat, n.years = dim(dat)[1])

    if (mode == 'simple')
    {
      TheLike <- as.numeric(all(whale.proj$pop$numbers > 1))
    }
    if (mode == 'fit')
    {
      TheLike <- exp(-1*(whale.likelihood(whale.proj))  + MLE)
    }

    Cumu <- Cumu + TheLike

    AveLike <- AveLike + TheLike

    Ntest <- Ntest +1
    while (Cumu > Threshold & Ndone < Nout) #check and see if cumulative likelihood is over threshold
    {
      show('wooo')

      Ndone <- Ndone + 1

      if (progbar == T)
      {
        setTxtProgressBar(pb, Ndone)
      }

      Cumu <- Cumu - Threshold
      Vals[Ndone,] <- data.frame(s.0.guess,s.rest.guess,K.guess,f.max.guess,
                                 last(whale.proj$pop$numbers),TheLike)
    }
  }

  Vals$AveLike <- AveLike/Ntest
  return(Vals)
}