arctan_conversion  <- function(x,lower_bound = 0, upper_bound = 1,invert = F)
{
  if (invert == T)
  {
    y <- lower_bound + (upper_bound - lower_bound)* ((atan(x) + pi/2)/pi)

  }
  if (invert == F)
  {
    y <-  tan((pi* (x - lower_bound)/(upper_bound - lower_bound)) - pi/2 )

  }

  return(y)

}