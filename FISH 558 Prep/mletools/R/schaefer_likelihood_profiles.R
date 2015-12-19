schaefer_profile <- function(r_range,k_range,msy_range,dat,k_guess = 3000, r_guess = .2,q_guess = .0003,sdev_guess = 0.1,error_type = 'observation',
                             analytic_q = T)
{

r_profile <- list()

for (r in seq_len(length(r_range)))
{

  mod_fit <- nlminb(start = c(k_guess,sdev_guess,q_guess), objective = schaefer_mle,dat = dat,
                                              use = 'optimization', pars_to_opt = c('k','sdev','q'), lower = c(0.01,.01,1e-7),
                                              upper = c(Inf,Inf,1), fixed_r = r_range[r],analytic_q = analytic_q,error_type = error_type)

  r_profile[[r]] <- data.frame(var = 'r',schaefer_mle(params = mod_fit$par,dat = dat,
                    pars_to_opt = c('k','sdev','q'), fixed_r = r_range[r],analytic_q = analytic_q, error_type = error_type)$model_summary)
}

  r_profile <- ldply(r_profile)


  k_profile <- list()
  for (i in seq_len(length(k_range)))
  {

    mod_fit <- nlminb(start = c(r_guess,sdev_guess,q_guess), objective = schaefer_mle,dat = dat,
                      use = 'optimization', pars_to_opt = c('r','sdev','q'), lower = c(0.01,.01,1e-7),
                      upper = c(1,Inf,1), fixed_k = k_range[i],analytic_q = analytic_q,error_type = error_type)

    k_profile[[i]] <- data.frame(var = 'k',schaefer_mle(params = mod_fit$par,dat = dat,
                                                        pars_to_opt = c('r','sdev','q'), fixed_k = k_range[i],analytic_q = analytic_q,
                                                        error_type = error_type)$model_summary)
  }

  k_profile <- ldply(k_profile)

  msy_profile <- list()
  for (m in seq_len(length(msy_range)))
  {
    mod_fit <- nlminb(start = c(r_guess,sdev_guess,q_guess), objective = schaefer_mle,dat = dat,
                      use = 'optimization', pars_to_opt = c('r','sdev','q'), lower = c(0.01,.01,1e-7),
                      upper = c(1,Inf,1), fixed_k = k_guess,analytic_q = analytic_q, fixed_MSY = msy_range[m],error_type = error_type)

    msy_profile[[m]] <- data.frame(var = 'msy',schaefer_mle(params = mod_fit$par,dat = dat,
                                                        pars_to_opt = c('r','sdev','q'), fixed_k = k_guess,analytic_q = analytic_q,fixed_MSY = msy_range[m],
                                                        error_type = error_type)$model_summary)
  }

  msy_profile <- ldply(msy_profile)



  r_profile$chi_square <- pchisq( 2* (r_profile$negative_log_likelihood - min(r_profile$negative_log_likelihood) ), df = 1)

  r_likelihood_plot <- (ggplot(r_profile, aes(r,negative_log_likelihood)) +
                          geom_point() +
                          geom_line(color = 'blue',alpha = 0.6) +
                          geom_hline(aes(yintercept = min(negative_log_likelihood) + 1.92)))

  r_chisqaure_plot <- (ggplot(r_profile, aes(r,chi_square)) +
                         geom_point() +
                         geom_line(color = 'blue',alpha = 0.6) +
                         geom_hline(yintercept = 0.95))




  k_profile$chi_square <- pchisq( 2* (k_profile$negative_log_likelihood - min(k_profile$negative_log_likelihood) ), df = 1)

  k_likelihood_plot <- (ggplot(k_profile, aes(k,negative_log_likelihood)) +
    geom_point() +
      geom_line(color = 'blue',alpha = 0.6) +
      geom_hline(aes(yintercept = min(negative_log_likelihood) + 1.92)))

  k_chisqaure_plot <- (ggplot(k_profile, aes(k,chi_square)) +
    geom_point() +
      geom_line(color = 'blue',alpha = 0.6) +
      geom_hline(yintercept = 0.95))

  msy_profile$chi_square <- pchisq( 2* (msy_profile$negative_log_likelihood - min(msy_profile$negative_log_likelihood) ), df = 1)

  msy_likelihood_plot <- (ggplot(msy_profile, aes(MSY,negative_log_likelihood)) +
                          geom_point() +
#                           geom_smooth(se = F) +
                            geom_line(color = 'blue',alpha = 0.6) +

                          geom_hline(aes(yintercept = min(negative_log_likelihood) + 1.92)))

  msy_chisqaure_plot <- (ggplot(msy_profile, aes(MSY,chi_square)) +
                         geom_point() +
#                          geom_smooth(se = F) +
                           geom_line(color = 'blue',alpha = 0.6) +

                         geom_hline(yintercept = 0.95))

  plots <- ls(pattern = '_plot')

  likelihood_profiles <- list()

  for (f in 1:length(plots))
  {
    eval(parse(text = paste('likelihood_profiles$',plots[f],'=',plots[f], sep = '')))
  }

return(likelihood_profiles)
}