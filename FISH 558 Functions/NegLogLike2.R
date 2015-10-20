NegLogLike2 <- function(Xcurr,dat,run_time)
{

#   a <- proc.time()
  out <- dat$blankpop

  out$s0 <- Xcurr$guess[Xcurr$var == 's0']

  out$q <- Xcurr$guess[Xcurr$var == 'q']

  out$lambda <- Xcurr$guess[Xcurr$var == 'lambda']

  egg.guess <- Xcurr$guess[grepl('egg',Xcurr$var)]

  out$eggs[2:(length(egg.guess)+1)] <- exp(egg.guess)

#   out$effort[1:length(dat$Effort)] <- dat$Effort
#
#   out$true.eggs[1:(length(dat$EggsObs))] <- dat$EggsObs
#
#   out$true.adults[1:(length(dat$adult.counts))] <- dat$adult.counts
#
#   out$true.interact[1:(length(dat$interaction.counts))] <- dat$interaction.counts


  for (t in 2:dat$run_time)
  {

    out$yearlings[t] <- out$eggs[t-1] * out$s0[1]

    out$juv[t] <- out$yearlings[t-1] * out$s.i[1] + (1-out$gamma[1])*out$juv[t-1]*out$s.j[1]

    out$adults[t] <- out$gamma[1]*out$juv[t-1]*out$s.j[1] +  out$adults[t-1] * out$s.a[1] -
      (1-out$lambda[1])*out$q[1]*out$effort[t-1]*out$adults[t-1]

  }

  egg_like <- sum(-dlnorm((out$true.eggs),meanlog = log(out$eggs), sdlog = dat$SDEgg, log = T), na.rm = T)

  agelike <- matrix(NA,nrow = 3,ncol = 1)

  for (t in 1:length(dat$n.years.data))
  {

    agelike[t] <- -dmultinom((dat$age.data[t,]),prob = as.numeric(out[dat$n.years.data[t],stages]), log = T)
  }

  age_structure_likes <- sum(agelike, na.rm = T)

  adult_like <- sum(-dlnorm((out$true.adults),meanlog = log(out$adults), sdlog = dat$log.sd.adults, log = T), na.rm = T)

  interact_like <- sum(-dlnorm(out$true.interact,meanlog = log(with(out,q*adults*effort)), sdlog = dat$log.sd.interaction, log = T), na.rm = T)

  total.like <- age_structure_likes + adult_like + interact_like + egg_like

  total.like[is.nan(total.like)] <- 10000

#   if (is.nan(total.like))
#   {
#     total.like <- 10000
#   }
#   proc.time() -a

  return(total.like)
}