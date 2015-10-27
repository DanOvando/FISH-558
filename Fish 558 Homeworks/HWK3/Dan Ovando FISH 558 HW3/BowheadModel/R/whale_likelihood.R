#' Calculate negative log likelhood of whale population
#'
#' \code{whale.likelihood} calculates negative log likelihood
#' @param whales an object created by make.whales()
#' @export
whale.likelihood <- function(whales)
{
    likelihood <- whales$pop %>%
    dplyr::select(-age.0) %>%
    mutate(predicted.abundance = rowSums(dplyr::select(.,contains('age.')))) %>%
    join(dplyr::select(whales$dat,year,abundance, abundance.sd,abundance.cv), by = 'year') %>%
    summarize(nll = sum( (log(predicted.abundance) - log(abundance))^2 / (2*(abundance.cv)^2), na.rm = T),
    l = prod(exp(-(log(predicted.abundance) - log(abundance))^2 / (2*log(abundance.cv)^2)), na.rm = T),
    nll2 = sum(-dnorm(log(abundance), log(predicted.abundance), sd = log(abundance.sd), log = T), na.rm = T))

  return(likelihood$nll)

}