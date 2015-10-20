#' Do a simple SIR to test validity of parameters
#'
#' \code{whale.optim} returns a posterior of calf and adult
#' survival, max fecundity, and K, such that the population
#' does not crash given an observed catch history
#' @param dat the raw whale data
#' @param s.0.
#' @param s.rest adult survival
#' @param K carrying capacity
#' @param f.max max fecundity
#' @export

whale.optim <- function(pars,dat,catch.dat)
{

  s.0 <- pars[1]

  s.rest <- pars[2]

  K <- pars[3]

  f.max <- pars[4]

  whales <- make.whales(dat = dat,catch.dat = catch.dat,f.max = f.max,s.0 = s.0,
                        s.rest = s.rest, extra.time = 0,Kinit = K)
  whale.proj <- whale.pop.model(whales = whales, dat = dat, n.years = dim(dat)[1])

  nll <- whale.likelihood(whale.proj)

  return(nll)

}