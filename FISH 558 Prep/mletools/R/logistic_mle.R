logistic_mle <- function(params,dat,index,n1,constrain = F,use = 'simulation')
{

  r <- params[1]

  k <- params[2]

  if (constrain == T)
  {
    r <- (atan(params[1]) + pi/2)/pi

    k <- n1 + (100000*n1 - n1)* ((atan(params[2]) + pi/2)/pi)
  }
  time <- dim(dat)[1]

  dat$observed <- dat[,index]

  dat$sdev <- dat[,paste('sd',index,sep = '_')]

  dat$predicted <- NA

  dat$predicted[1] <- n1

  for (t in 2:time)
  {

    pastpop <- dat$predicted[t-1]

    dat$predicted[t] <- pmin(k,pmax(0,pastpop + pastpop*r*(1-pastpop/k)))

  }
  dat$nll <- with(dat,log(sdev) + 0.5 * log(2 * pi) + (observed - predicted) ^ 2 / (2 * sdev ^ 2))

  dat$r <- r

  dat$k <- k

  nll <- sum(dat$nll, na.rm = T)

  if (use == 'optimization')
  {
    output <- nll
  }

  if (use == 'simulation')
  {
    output <- dat
  }

  return(output)
}