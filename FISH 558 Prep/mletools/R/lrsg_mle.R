lrsg_mle <- function(params, fixed_b0 = NA, fixed_s = NA,
                     fixed_sdev = NA, fixed_rlags = 4,fixed_z = NA, dat,constrain = F,use = 'simulation',
                     pars_to_opt = c('b0','s','z','sdev','q'),analytic_q = F)

{
  for (p in 1:length(pars_to_opt))
  {
    eval(parse( text = paste(pars_to_opt[p],'<-',params[p], sep = '')))
  }

  not_opted <- c('b0','s','z','sdev','q','rlags')[!(c('b0','s','z','sdev','q','rlags') %in% pars_to_opt)]

  for (n in seq_len(length(not_opted)))
  {
    eval(parse( text = paste(not_opted[n],'<- fixed_',not_opted[n], sep = '')))
  }


  dat <- lrsg_model(dat = dat, b0 = b0, s = s, z = z, q = q, rlags = rlags)$dat

  dat$nll <- with(dat,log(sdev) + 0.5 * log(2 * pi) + (log(predicted_cpue) - log(cpue)) ^ 2 / (2 * sdev ^ 2))

  dat$s <- s

  dat$z <- z

  dat$q <- q

  dat$rlags <- rlags

  dat$b0 <- b0

  dat$sdev <- sdev

  nll <- sum(dat$nll, na.rm = T)

  model_summary <- data.frame(negative_log_likelihood = nll,
                              likelihood = -exp(nll),
                              b0 = b0,
                              z = z,
                              s = s,
                              sdev = sdev,
                              rlags = rlags,
                              q = q,
                              MSY = unique(dat$MSY))

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