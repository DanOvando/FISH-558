#' Run a SIR algorithm on the supplied whale population
#'
#' \code{whale.sir} returns a posterior of calf and adult
#' survival, max fecundity, and K, given the supplied likelihood function
#' @param Nout the number of desired draws from the posterior
#' @param dat the raw whale data
#' @param s.0.lower lower bound of calf survival
#' @param s.0.upper upper bound of calf survivaletc
#' with the rest of the lower and uppers
#' @param mode 'simple' or 'fit' tells it to either run
#' SIR on crashing th population, 'fit' tells it to fit the
#' model
#' @param use.catch T or F to add catch to model
#' @param MLE the guess of the minimum of the NLL from MLE
#' @param mode tells the model whether to fit binary SIR or
#' full likelihood SIR
#' @export

whale.sir <- function(Nout=1000,dat,catch.dat = catch.dat,s.0.lower = 0.8,s.0.upper = 1,
                      s.rest.lower = 0.9, s.rest.upper = 1,
                      k.lower = 10000,k.upper = 20000,
                      f.max.lower = 0.25,f.max.upper = 0.33,progbar = F,mode = 'simple', MLE = 2.22
                      , use.catch = T)
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

  Threshold <- exp(-MLE)
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
                          s.rest = s.rest.guess, K = K.guess, use.catch = use.catch)

    #Run the model using the created object
    whale.proj <- whale.pop.model(whales = whales, dat = dat, n.years = dim(dat)[1])

    if (mode == 'simple')
    {
      TheLike <- as.numeric(all(whale.proj$pop$numbers > 1)) #binary likelihood
    }
    if (mode == 'fit')
    {
#       TheLike <- exp(-1*(whale.likelihood(whale.proj) - MLE)) # scale likelihood by the estimate of best likelihood
      TheLike <- exp(-1*(whale.likelihood(whale.proj))) #calculate NLL and convert to likelihood

    }

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
      Vals[Ndone,] <- data.frame(s.0.guess,s.rest.guess,K.guess,f.max.guess,
                                 last(whale.proj$pop$numbers),TheLike)
    }

  }

  Vals$AveLike <- AveLike/Ntest
  return(Vals)
}