#' Convert moment to shape parameters
#'
#' \code{moments_to_shape} takes the mean and cv
#' of a beta distribution and converts it to the shape
#' parameters
#' @param mu the mean of the beta distribution
#' @param cv the cv of the beta distribution

moments_to_shape <- function(mu,cv)
{

  var <- (mu*cv)^2 #calculate variance

  a <- ( (1-mu)/var - 1/mu )*mu^2

  b <- a*(1/mu - 1)

  if (mu >=1){warning('mu is greater than or equal to 1, reduce it!')}

  return(list(shape1 = a, shape2 = b))
}
