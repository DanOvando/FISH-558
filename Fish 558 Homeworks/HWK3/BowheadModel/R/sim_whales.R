#' Simulate whale populations
#'
#' \code{sim.whales} projects the whale population forward
#' using outputs from the SIR algorithm (or for any matrix of
#' possible params)
#' @param  dat the raw data
#' @param catch.dat the raw catch data
#' @param possibe.whales a data.frame of possible parameters
#' for whale model
#' @param extra.time the number of years to pass forward
#' @export

sim.whales <- function(dat,catch.dat,possible.whales,extra.time = NA, extra.catch = NA, use.catch = T)
{

  simmed.whales <-  foreach(i=1:dim(possible.whales)[1], .combine = rbind) %do%
  {
    pars <- possible.whales[i,]

    whale.fit <- make.whales(dat = dat, catch.dat = catch.dat, s.0 = pars$s.0,
                             s.rest = pars$s.rest, f.max =pars$f.max,
                             K = pars$K, extra.time = extra.time, extra.catch = extra.catch,
                             use.catch = use.catch)

    fitted.whales <- whale.pop.model(whales = whale.fit, dat = dat, n.years = dim(whale.fit$pop)[1])

    whale.trajectory <- (fitted.whales$pop) %>%
      select(-age.0, -numbers) %>%
      gather('age','numbers',age.1:age.13, convert = T) %>%
      group_by(year) %>%
      summarize(predicted.whales = sum(numbers, na.rm = T),
                predicted.mature.whales = numbers[age == paste('age',whales$life$max.age, sep = '.')],
                predicted.calves = numbers[age == 'age0']) %>%
      left_join(dat, by = 'year') %>%
      select(-catch) %>%
      left_join(whale.fit$catch.dat, by = 'year')

    whale.run <- suppressWarnings(data.frame(whale.trajectory, pars,iteration = i, catchlevel = extra.catch))

    return(whale.run)
  }

  return(simmed.whales)
}