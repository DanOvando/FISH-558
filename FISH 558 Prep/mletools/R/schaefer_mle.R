schaefer_mle <- function(params,fixed_r = NA,fixed_k = NA, fixed_q = NA,
                         fixed_sdev = NA, dat,constrain = F,use = 'simulation',
                         pars_to_opt = c('r','k','q','sdev'),analytic_q = F, fixed_MSY = NA,
                         error_type = 'observation')
{
  for (p in 1:length(pars_to_opt))
  {
    eval(parse( text = paste(pars_to_opt[p],'<-',params[p], sep = '')))
  }

  not_opted <- c('r','k','q','sdev')[!(c('r','k','q','sdev') %in% pars_to_opt)]

  for (n in seq_len(length(not_opted)))
  {
    eval(parse( text = paste(not_opted[n],'<- fixed_',not_opted[n], sep = '')))
  }

  if (constrain == T)
  {
    r <- (atan(params[1]) + pi/2)/pi

    k <- n1 + (100000*n1 - n1)* ((atan(params[2]) + pi/2)/pi)
  }

  if (is.na(fixed_MSY) == F)
  {
    k <- (4*fixed_MSY)/r
  }


  time <- dim(dat)[1]

  dat$predicted_n <- NA

  dat$predicted_n[1] <- k

  if (error_type == 'process')
  {
    dat$predicted_n[1] <- dat$cpue[1]/q
  }

  dat$predicted_cpue <- NA

  dat$predicted_cpue[1] <- q *  dat$predicted_n[1]

  for (t in 2:time)
  {

    pastpop <- dat$predicted_n[t-1]

    if (error_type == 'observation')
    {

    dat$predicted_n[t] <- pmin(k,pmax(0.001,pastpop + pastpop*r*(1-pastpop/k) - dat$catch[t-1]))

    dat$predicted_cpue[t] <- q * dat$predicted_n[t]
    }
    if (error_type == 'process')
    {
      index_pop <- 1/q * dat$cpue[t-1]

      dat$predicted_n[t] <- pmin(k,pmax(0.001,index_pop  + index_pop*r*(1- 1/k * index_pop) - dat$catch[t-1]))

      dat$predicted_cpue[t] <- q * dat$predicted_n[t]

    }
  }

  if (analytic_q == T)
  {
    q <- exp(1/time * sum(log(dat$cpue) - log(dat$predicted_n)  ))

    dat$predicted_cpue <- q * dat$predicted_n
  }

  dat$nll <- with(dat,log(sdev) + 0.5 * log(2 * pi) + (log(predicted_cpue) - log(cpue)) ^ 2 / (2 * sdev ^ 2))

  MSY <- (r*k)/4

  dat$r <- r

  dat$k <- k

  dat$MSY <- MSY

  dat$q <- q

  dat$sdev <- sdev

  nll <- sum(dat$nll, na.rm = T)


  model_summary <- data.frame(negative_log_likelihood = nll,
    likelihood = exp(-nll),
     r = r,
    k = k,
    q = q,
    sdev = sdev,
    MSY = MSY)

  if (use == 'optimization')
  {
    output <- nll
  }

  if (use == 'simulation')
  {
    output <- list(dat = dat, model_summary = (model_summary))
  }

  return(output)
}