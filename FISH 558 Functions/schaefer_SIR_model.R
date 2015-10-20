schaefer.model <- function(pars,dat,time,mode = 'optim')
{
  #   r,k,phi

  r <- exp(pars[1])

  k <- exp(pars[2])

  phi <- exp(pars[3])

  time <- dim(dat)[1]

  dat$k <- k #runif(1,500,8000)

  dat$phi <- phi #runif(1,0.8,1.2)

  dat$r <- r #runif(1,.15,.45)

  dat$biomass[1] <- dat$phi[1] * dat$k[1]

  for (t in 2:time)
  {

    dat$biomass[t] <- pmax(1e-8,with(dat,biomass[t-1] + biomass[t-1] * r[t-1] * (1 - biomass[t-1]/k[t-1]) - catch[t-1]))

  }

  dat$q <- exp(1/time * sum(  log(dat$catch.rate / dat$biomass)))

  dat$index <- dat$biomass * dat$q

  dat$sigma <- sqrt(1/time * sum( (log(dat$catch.rate) - log(dat$q*dat$biomass))^2))

  dat$thelike <- with(dat,{
    (1/(sqrt(2*pi)*sigma*catch.rate)* exp(- (log(catch.rate) - log(index))^2 / (2*sigma^2)))
    #     log(sigma * catch.rate) + (log(catch.rate) - log(index))^2 / (2*sigma^2)

  }
  )

  dat$total_like <- prod(dat$thelike)

  output <- dat
  if (mode == 'optim')
  {
    output <- as.numeric(dat$total_like[1])
  }

  return(output)
}