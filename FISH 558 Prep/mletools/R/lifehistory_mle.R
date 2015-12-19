lifehistory_mle <- function(params,dat,n1,index,constrain = F,use = 'simulation'
                            , Ntarget = 2,gamma = 2,M = 100,target_N = F,rainfall = 150,fityear = 1978)

{

  a <- params[1]

  b <- params[2]

  g <- params[3]

  f <- params[4]

  if (constrain == T)
  {
    r <- (atan(params[1]) + pi/2)/pi

    k <- n1 + (100000*n1 - n1)* ((atan(params[2]) + pi/2)/pi)
  }

  time <- dim(dat)[1]

  dat$total_food <- NA

  dat$food_per_animal <- NA

  dat$births <- NA

  dat$predicted_calf_survival <- NA

  dat$predicted_adult_survival <- NA

  dat$n_predicted <- NA

  dat$n_predicted2 <- NA

  dat$n_predicted2[1] <- n1

  dat$n_predicted[1] <- n1

  dat$observed <- dat[,index]

  dat$sdev <- dat[,paste('sd',index,sep = '_')]

  for (t in 2:(time+1))
  {

    pastpop <- dat$n_predicted[t-1]

    dat$total_food[t-1] <- 1.25* dat$dry_season_rain[t-1]

    dat$food_per_animal[t-1] <-  (dat$total_food[t-1] * 1000000 )/pastpop

    dat$births[t-1] <- 0.4*pastpop

    dat$predicted_calf_survival[t-1] <- (a*dat$food_per_animal[t-1])/(b + dat$food_per_animal[t-1])

    dat$predicted_adult_survival[t-1] <- (g*dat$food_per_animal[t-1])/(f + dat$food_per_animal[t-1])

    if (t <= time)
    {
            dat$n_predicted[t] <- pmax(0,pastpop*dat$predicted_adult_survival[t-1] + dat$births[t-1]*dat$predicted_calf_survival[t-1])

#       dat$n_predicted[t] <- pmax(0,pastpop * ( (g * 1000000*1.25* dat$dry_season_rain[t-1] / pastpop) / (f + 1000000*1.25* dat$dry_season_rain[t-1] / pastpop)) +
#                                    0.4*pastpop * ( (a * 1000000*1.25* dat$dry_season_rain[t-1] / pastpop) / (b + 1000000*1.25* dat$dry_season_rain[t-1] / pastpop)))
#
       }
  }

  dat$predicted_monthly_mortality <- .25*(1-(dat$predicted_adult_survival))

  dat$nll_census <- with(dat,((observed) - (n_predicted)) ^ 2 / (2 * sdev ^ 2))

  dat$nll_calf_survival <- with(dat,((calf_survival) - (predicted_calf_survival)) ^ 2 / (2 * sd_calf_survival ^ 2))

  dat$nll_adult_mortality <- with(dat,((adult_mort) - (predicted_monthly_mortality) ) ^ 2 / (2 * sd_adult_mort ^ 2))

  dat$a <- a

  dat$b <- b

  dat$g <- g

  dat$f <- f

  nll_census <- sum(dat$nll_census[dat$year <= fityear], na.rm = T)

  nll_calf_survival <- sum(dat$nll_calf_survival[dat$year <= fityear], na.rm = T)

  nll_adult_mortality <- sum(dat$nll_adult_mortality[dat$year <= fityear], na.rm = T)

  a_prime <- b*f

  b_prime <- (1.25*rainfall * 1000000) * (b + f - g*b - 0.4*a*f)

  c_prime <- (1.25*rainfall * 1000000)^2 * (1 - g - 0.4*a)

  n_eq <- (-b_prime + sqrt((b_prime)^2 - 4 * a_prime * c_prime))/(2*a_prime)

  dat$n_eq <- last(dat$n_predicted)

  dat$eq_rain <- rainfall

  if (use == 'optimization')
  {
    output <- sum(nll_census,nll_calf_survival,nll_adult_mortality)

    if (target_N == T)
    {
      output <- output + ((last(dat$n_predicted) - Ntarget * n1)^gamma / M )
#       output <- output + (n_eq - Ntarget * n1)^gamma / M

    }
  }

  if (use == 'simulation')
  {
    output <- dat
  }

  return(output)
}