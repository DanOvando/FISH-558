#' Process whale data
#'
#' \code{make.whales} returns a last object containing all
#' the things you might need for a whale function
#' @param dat the raw data
#' @param catch the catch data if any. can be NA to not add any catch
#' @param default.catch the amount of catch to put in in any year in which
#' no catch data exists
#' @param s.0 calf survival
#' @param s.rest adult survival
#' @param f.max max fecundity
#' @param f.0 fecundity at unfished EQ
#' @param z compensation ratio
#' @param Kinit population in age 1, used to scale carrying capacity
#' @param max.age plus group age
#' @param extra.time number indicating how many years beyond the last
#' year of data to add

make.whales <- function(dat,catch.dat,default.catch = 0,s.0 = 0.4, s.rest = 0.8, f.max =3, z = 2.39, Kinit = 10000, max.age = 13, extra.time = 0)
{

  K <- Kinit * exp((0:(max.age - 2))* log(s.rest))

  K <- c(K,  (last(K) * exp(log(s.rest))) / (1-exp(log(s.rest))))   # Determine numbers at age

  f.0 <- (sum(K) * (1-s.rest)) / (last(K)*s.0) #find f.0, based on births = deaths at EQ

  whales <- list()

  whales$dat <- dat
  #store whale things
  whales$life <- list(s.0 = s.0, s.rest = s.rest, f.0 = f.0, f.max = f.max, z = z, K = K, max.age = max.age)

  whales$dat$abundance.sd <- whales$dat$abundance * whales$dat$abundance.cv

  if (extra.time >0) #tack on extra years if needed
  {
    time <- c(dat$year, last(dat$year) + 1:extra.time)
  } else{
    time <- dat$year
  }

  # Convert age and catch to matrix format for faster processing
  age.classes <- as.data.frame(matrix(NA, nrow = length(time), ncol = max.age + 1))

  colnames(age.classes) <- paste('age', 0:max.age, sep = '.')

  age.classes[1,2:dim(age.classes)[2]] <- K

  age.classes[1,1] <-  last(K) * f.0  # add in starting conditions calves

  catch.at.age <- as.data.frame(matrix(default.catch, nrow = length(time), ncol = max.age + 1)) #catch at age can be set to custom value

  colnames(catch.at.age) <- paste('catch', 0:max.age, sep = '.')

  if (class(catch.dat) != 'logical') #add in observed catch data if needed
  {
    catch.at.age$year <- time

    catch.at.age<- left_join(catch.at.age,catch.dat, by = 'year')

    catch.at.age[, paste('catch', 1:max.age, sep = '.')] <- catch.at.age$catch / length(1:max.age)

    catch.at.age <- dplyr::select(catch.at.age, -year, -catch)
  }

  whales$pop <- data.frame(year = time, numbers = rep(NA,length(time)),age.classes, catch.at.age) #make population object

  # A whole bunch of stuff to store indexing
  whales$adults <- which(grepl('age.', colnames(whales$pop)) & colnames(whales$pop) != 'age.0')

  whales$mid.adults <-  which(grepl('age.', colnames(whales$pop)) & colnames(whales$pop) != 'age.0' & colnames(whales$pop) != 'age.1'
                              & colnames(whales$pop) != paste('age',max.age, sep = '.'))

  whales$plus.adults <- which(colnames(whales$pop) == paste('age',max.age, sep = '.'))

  whales$pop <- left_join(whales$pop,dat[,c('year','catch')], by = 'year')

  whales$pop[, paste('catch', 1:max.age, sep = '.')] <- whales$pop$catch / length(1:max.age)

  whales$catch.loc <- paste('catch',1:(whales$life$max.age-2), sep = '.')

  whales$next.to.last <- paste('age',whales$life$max.age-1, sep = '.')

  whales$next.to.last.catch <- paste('catch',whales$life$max.age-1, sep = '.')

  whales$max.age <- paste('age',whales$life$max.age, sep = '.')

  whales$max.catch <- paste('catch',whales$life$max.age, sep = '.')

  whales$age.matrix <- as.matrix(age.classes)

  whales$catch.matrix <- as.matrix(catch.at.age)

  # Idea for gathering on multiple groups of data
  #   whales$pop %>%
  #     gather(key, value, -id, -time) %>%
  #     extract(key, c("question", "loop_number"), "(Q.\\..)\\.(.)") %>%
  #     spread(question, value)

  return(whales)
}

#' Simulate whale population
#'
#' \code{whale.pop.model} Projects a whale population forward using input from
#' make.whales
#' @param whales an object created by make.whales()
#' @param dat the raw data
#' @param nyears the number of years to simulate

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

    pop[t, mid.adults] <- (pop[t-1,mid.adults-1] - whales$catch.matrix[t-1,mid.adults - 1]) *
      whales$life$s.rest #non plus group

    pop[t,max.age] <- ((pop[t-1,max.age - 1] -
                          whales$catch.matrix[t-1,max.age - 1]) *
                         whales$life$s.rest +
                         (pop[t-1,max.age] -
                            whales$catch.matrix[t-1,max.age]) *
                         whales$life$s.rest ) #plus group

    pop[t,1] <- last(pop[t,]) * (whales$life$f.0 + (whales$life$f.max - whales$life$f.0) *
                                   (1 - ( sum(pop[t,adults]) / sum(whales$life$K) )^whales$life$z)) #make calves

  }

  #   whales$pop <- pop
  whales$pop[,grepl('age.',colnames(whales$pop))] <- as.data.frame(pop)
  return(whales)
}

whale.likelihood <- function(whales)
{

  likelihood <- whales$pop %>%
    mutate(predicted.abundance = rowSums(dplyr::select(.,contains('age.')))) %>%
    join(dplyr::select(whales$dat,year,abundance, abundance.sd,abundance.cv), by = 'year') %>%
    summarize(nll = sum( (log(predicted.abundance) - log(abundance))^2 / (2*log(abundance.cv)^2), na.rm = T))

  return(likelihood$nll)

}




