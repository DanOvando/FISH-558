tune_scalar <- function(scalar,n_sim = 10000,par_init,parnames,dat,vcov,n_burn=0,n_thin=1,
                     vcov_augment = 1,prog_bar = T,jumpyness = 1, seed = NA,targ_accept_rate = 0.3)
{

  mcmc_results <- mozzy_mcmc(par_init = par_init,
                             parnames = parnames, dat = dat,
                             vcov = vcov, prog_bar = prog_bar, n_sim = n_sim, n_burn = round(.5*n_sim,0),
                             n_thin = 1,vcov_augment = scalar,jumpyness = .001,targ_accept_rate = targ_accept_rate)

  return((mcmc_results$accepted_runs$perc_selected - targ_accept_rate)^2)

}