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
#' @export

make.whales <- function(dat,catch.dat,default.catch = 0,s.0 = 0.4, s.rest = 0.8, f.max =3, z = 2.39, K = 10000, max.age = 13, extra.time = 0)
{

  #   K <- Kinit * exp((0:(max.age - 2))* log(s.rest))
  #
  #   K <- c(K,  (last(K) * exp(log(s.rest))) / (1-exp(log(s.rest))))   # Determine numbers at age
  init.calfs <- K / sum(c(s.0 * s.rest^(0:11), (s.0*s.rest^12) / (1 - s.rest)))

  init.adults <- init.calfs * c(s.0 * s.rest^(0:11), (s.0*s.rest^12) / (1 - s.rest))

  f.0 <- ((K) * (1-s.rest)) / (last(init.adults) * s.0) #find f.0, based on births = deaths at EQ

#   f.0 <- (1-s.rest)/(s.0*s.rest^12)

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

  age.classes[1,2:dim(age.classes)[2]] <- init.adults

  age.classes[1,1] <-  init.calfs  # add in starting conditions calves

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