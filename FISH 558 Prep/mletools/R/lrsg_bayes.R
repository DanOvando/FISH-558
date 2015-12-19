lrsg_bayes <- function(i = 0,b0_range, s_range, z_range,rlags = 4,dat)
{
  b0 <- runif(1,b0_range[1],b0_range[2])

  s <- runif(1,s_range[1],s_range[2])

  z <- runif(1,z_range[1],z_range[2])

  time <- dim(dat)[1]

  r0 <- b0 * (1-s)

  a <- b0/r0 * (1 - (z - 0.2)/(0.8*z))

  b <- (z- 0.2)/(0.8 * z * r0)

  dat$predicted_b[1] <- b0

  for (t in 1:time)
  {
    if ((t - rlags)<=0){
      dat$predicted_recruitment[t] <- r0
    }else {
      dat$predicted_recruitment[t] <- dat$predicted_b[t-rlags] / (a +  b*dat$predicted_b[t-rlags])
    }
    if (t < time){
      dat$predicted_b[t+1] <- pmax(.01,s*dat$predicted_b[t] + dat$predicted_recruitment[t] - dat$catch[t])
    }

  }

  #   dat$lag_predicted_b <- lag(dat$predicted_b,rlags)

  dat$bmsy <- 1/b * sqrt( a / (1-s) - a)

  dat$MSY <- dat$bmsy * ( s - 1 + (1 / (a + b*dat$bmsy)))
if (any(dat$MSY < 0)){ browser()}
  opt_q_sdev <- function(params,dat){

    sdev <- params[1]

    q <- params[2]

    predicted_cpue <- dat$predicted_b * q

    nll <- log(sdev) + 0.5 * log(2 * pi) + (log(predicted_cpue) - log(dat$cpue)) ^ 2 / (2 * sdev ^ 2)

    return(sum(nll, na.rm = T))
  }

  q_and_sdev <- nlminb(start = c(sdev_guess, q_guess), objective = opt_q_sdev, dat = dat,
                       lower = c(1e-7,1e-7), upper = c(Inf,1))

  dat$sdev <-  q_and_sdev$par[1]

  dat$q <- q_and_sdev$par[2]

  dat$predicted_cpue <- dat$predicted_b * dat$q

  dat$nll <- with(dat,log(sdev) + 0.5 * log(2 * pi) + (log(predicted_cpue) - log(cpue)) ^ 2 / (2 * sdev ^ 2))

  total_nll <- sum(dat$nll, na.rm = T)

  it_summary <- data.frame(
    iteration = i,
    negative_log_likelihood = total_nll,
    likelihood = exp(-total_nll),
    b0 = b0,
    bmsy = mean(dat$bmsy),
    MSY = mean(dat$MSY),
    s = s,
    z = z,
    rlags = rlags,
    q = mean(dat$q),
    sdev = mean(dat$sdev)
  )

  return(it_summary)
}