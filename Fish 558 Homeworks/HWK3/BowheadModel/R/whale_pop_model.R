#' Simulate whale population
#'
#' \code{whale.pop.model} Projects a whale population forward using input from
#' make.whales
#' @param whales an object created by make.whales()
#' @param dat the raw data
#' @param nyears the number of years to simulate
#' @export

whale.pop.model <- function(whales,dat,n.years)
{

  pop <- whales$age.matrix

  max.age <- dim(pop)[2]

  #indexing outside of loop
  adults <- 2 : max.age

  mid.adults <- 3 : (max.age-1)

  max.age <- dim(pop)[2]

  for (t in 2:n.years)
  {

    pop[t,2] <- pop[t-1,1] * whales$life$s.0 #calf survival

    pop[t, mid.adults] <- pmax(1e-3,(pop[t-1,mid.adults-1] - whales$catch.matrix[t-1,mid.adults - 1]) *
                                 whales$life$s.rest) #non plus group

    pop[t,max.age] <- pmax(1e-3,((pop[t-1,max.age - 1] -
                                    whales$catch.matrix[t-1,max.age - 1]) *
                                   whales$life$s.rest)) +
      pmax(1e-3,(pop[t-1,max.age] -
                   whales$catch.matrix[t-1,max.age]) *
             whales$life$s.rest) #plus group

    pop[t,1] <- pmax(1e-3,last(pop[t,]) * (whales$life$f.0 + (whales$life$f.max - whales$life$f.0) *
                                             (1 - ( sum(pop[t,adults]) / (whales$life$K) )^whales$life$z))) #make calves

  }
  whales$pop[,grepl('age.',colnames(whales$pop))] <- as.data.frame(pop)

  whales$pop$numbers <- rowSums(  pop[,grepl('age.',colnames(pop))])


  return(whales)
}
