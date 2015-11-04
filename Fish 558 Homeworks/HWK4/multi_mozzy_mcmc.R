#' Run multi-chain MCMC on mosquitoes
#'
#' \code{multi_mozzy_mcmc} runs potentially multiple chains
#'  MCMC in parallel
#' mosquito count data
#' @param par_guess a dataframe of starting parameter vectors
#' @param num_starts the number of different starting guesses you wan
#' @param jitfactor the multiplier for the jitter function
#' @param numcores number of cores to use in parallel
#' @param varnames names of the parameters
#' @param dat are the raw stream-count-distance data
#' @param covar is the variance-covariance matrix of pars
#' @param n_sim number of simulations including burn in
#' @param n_burn number of simulations to burn in
#' @param n_thin the number of iterations to thin by
#' @param vcov_augment factor to multiply vcov by in the burn in period
#' @param prog_bar T or F to output progress bar

multi_mozzy_mcmc<-function(par_guess,num_starts = 1, jitfactor = 1,numcores = 1,parnames,dat,vcov,n_sim=1000,n_burn=0,n_thin=1,
                           vcov_augment = 1,prog_bar = F,jumpyness = 1)
{

  registerDoMC(cores=numcores)

    jitter_start <- function(i,pars, jitfactor)
  {

    newpars <- pars

    newpars[1,] <- jitter(as.numeric(newpars), factor = jitfactor)

    return(newpars)
  }

  if (num_starts > 1) {
starting_pars <- lapply(1:num_starts, jitter_start, pars = par_guess,jitfactor = jitfactor) }else{
  starting_pars <- lapply(1, jitter_start, pars = par_guess, jitfactor = 0)
}


  multi_chain <- foreach(i = 1:length(starting_pars)) %dopar%
     mozzy_mcmc(par_init = starting_pars[[i]], parnames = colnames(par_guess), dat = dat,
                               vcov = vcov, prog_bar = prog_bar, n_sim = n_sim, n_burn = n_burn, n_thin = n_thin,
                               vcov_augment = vcov_augment,seed = runif(1,1,1000))


  posteriors <- lapply(seq(along = multi_chain), function(i)    data.frame(run = i, multi_chain[[i]]$posteriors)) %>%
    ldply()

return(list(posteriors = posteriors,multi_chain = multi_chain))

}