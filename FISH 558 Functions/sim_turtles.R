sim_turtles <- function(i,mc.params,dat,run_time)
{
  Xcurr <- data.frame((mc.params[i,]))

#   colnames(Xcurr) <- c('var','guess')

  out <- dat$blankpop

  out$s0 <- Xcurr$s0

  out$q <- Xcurr$q

  out$lambda <- Xcurr$lambda

  egg.guess <- Xcurr[,grepl('egg',colnames((Xcurr)))]

  out$eggs[2:(length(egg.guess)+1)] <- as.numeric(exp(egg.guess))

  out$effort[is.na(out$effort)] <- Xcurr$future.effort

  out$future.effort <- Xcurr$future.effort

  for (t in 2:run_time)
  {

    if (is.na(out$eggs[t]))
    {
    out$eggs[t] <- 15.5 * out$adults[t-1]
    }

    out$yearlings[t] <- out$eggs[t-1] * out$s0[1]

    out$juv[t] <- out$yearlings[t-1] * out$s.i[1] + (1-out$gamma[1])*out$juv[t-1]*out$s.j[1]

    out$adults[t] <- pmax(1e-5,out$gamma[1]*out$juv[t-1]*out$s.j[1] +  out$adults[t-1] * out$s.a[1] -
      (1-out$lambda[1])*out$q[1]*out$effort[t-1]*out$adults[t-1])

  }
  return(out)
}