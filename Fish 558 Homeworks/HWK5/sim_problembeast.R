#' runs problembeast model
#'
#'\code{sim_problembeast} simulates populations of
#'the endangered problembeast forward under
#'different states of nature and hunting
#'@param pb a list of problembeast population
#'parameters
#'@param sim_years number of years to run
#'@param hunted_patch the patch in which hunting
#'occurs
#'@param s_j juvenile survival
#'@param _value the price per hunted animal
#'at various ages
#'@param results blank matrix to store results
#'@param recruits index of recruits
#'@param juveniles index of juveniles
#'@param adults index of all adults
#'@param mature index of mature adults
#'@param plus index of the plus group
#'
sim_problembeast <- function(pb,sim_years = 101,hunted_patch = 0,s_j = .75,
                             recruit_value = 1, juv_value = 5, adult_value = 30,
                             results, recruits,juveniles,adults,mature,plus)
{

pb$s_j <- s_j

# resultnames <- c('year','adults','hunted','hunt_value','hunted_patch','s_j')
#
# results <- as.data.frame(matrix(0,nrow = sim_years+1,ncol = length(resultnames)))
#
# colnames(results) <- resultnames

results$year <- 1:(sim_years+1)

results$hunted_patch <- hunted_patch

results$s_j <- s_j

# recruits <- which(colnames(pb$n_at_age) %in% 'age.0')
#
# juveniles <- which(colnames(pb$n_at_age) %in% paste('age',1:8,sep = '.'))
#
# adults <- which(colnames(pb$n_at_age) %in% paste('age',9:14,sep = '.'))
#
# mature <- which(colnames(pb$n_at_age) %in% paste('age',9:15,sep = '.'))
#
# plus <- which(colnames(pb$n_at_age) == 'age.15')

phi <- matrix(1, nrow = 3, ncol = sim_years)

if (hunted_patch>0){

  shape1 <- eval(parse(text = paste('pb$phi_',hunted_patch,'_shape',sep = '') ))$shape1

  shape2 <- eval(parse(text = paste('pb$phi_',hunted_patch,'_shape',sep = '') ))$shape2

  phi[hunted_patch,] <-  rbeta(sim_years,shape1 = shape1, shape2 = shape2)

}

n_adults <- sum(pb$n_at_age[1, mature])

results$adults[1] <- n_adults

pb$n_at_age <- as.matrix(pb$n_at_age)

  for (y in 1:(sim_years))
  {

    hunted <- 0

    juv_survive <- rbinom(length(juveniles),(pb$n_at_age[y,juveniles - 1]),pb$s_j)

    juv_survive_hunting <- rbinom(length(juveniles),(juv_survive),phi[2,y])

    juv_hunted <- juv_survive - juv_survive_hunting

    pb$n_at_age[y + 1,juveniles] <- juv_survive_hunting

    adults_survive <- rbinom(length(adults),(pb$n_at_age[y,adults - 1]),pb$s_a)

    adults_survive_hunting <- rbinom(length(adults), (adults_survive),phi[3,y])

    adults_hunted <- adults_survive - adults_survive_hunting

    pb$n_at_age[y + 1, adults] <- adults_survive_hunting

    enter_plus_survive <- rbinom(length(plus),(pb$n_at_age[y,plus - 1]),pb$s_a)

    enter_plus_survive_hunting <- rbinom(length(plus), (enter_plus_survive),phi[3,y])

    enter_plus_hunted <- enter_plus_survive - enter_plus_survive_hunting

    enter_plus <- enter_plus_survive_hunting

    stay_plus_survive <- rbinom(length(plus),(pb$n_at_age[y,plus]),pb$s_a)

    stay_plus_survive_hunted <- rbinom(length(plus), (stay_plus_survive),phi[3,y])

    stay_plus_hunted <- stay_plus_survive - stay_plus_survive_hunted

    stay_plus <- stay_plus_survive_hunted

    hunted_in_two <- sum(juv_hunted)

    hunted_in_three <- sum(adults_hunted) + sum(enter_plus_hunted) + sum(stay_plus_hunted)

    hunted <-  hunted_in_two + hunted_in_three

    hunted_value <- hunted_in_two * juv_value + hunted_in_three * adult_value

    pb$n_at_age[y + 1, plus] <- enter_plus + stay_plus

    n_adults <- sum(pb$n_at_age[y + 1, mature])

    recruits_survive <- rbinom(length(recruits), n_adults, pb$b)

    recruits_survive_hunting <- rbinom(length(recruits), (recruits_survive),phi[1,y])

    recruits_hunted <- recruits_survive - recruits_survive_hunting

    if (hunted_patch >1){
    results$hunted[y] <- hunted
    results$hunt_value[y] <- hunted_value
    }
    if(hunted_patch == 1)
    {
      results$hunted[y + 1] <- recruits_hunted
      results$hunt_value[y+1] <- recruits_hunted * recruit_value
    }

    pb$n_at_age[y + 1,recruits] <- recruits_survive_hunting

    results$adults[y+1] <- n_adults

  }

results <- results[1:sim_years,]

return(results)

}