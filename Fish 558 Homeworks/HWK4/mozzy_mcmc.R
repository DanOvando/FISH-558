#' Run MCMC on mosquitoes
#'
#' \code{Mozzy_MCMC} runs an MCMC on
#' mosquito count data
#' @param par_init initival vector of parameter guesses
#' @param varnames names of the parameters
#' @param dat are the raw stream-count-distance data
#' @param covar is the variance-covariance matrix of pars
#' @param n_sim number of simulations including burn in
#' @param n_burn number of simulations to burn in
#' @param n_thin the number of iterations to thin by
#' @param vcov_augment factor to multiply vcov by in the burn in period
#' @param prog_bar T or F to output progress bar

mozzy_mcmc<-function(par_init,parnames,dat,vcov,n_sim=1000,n_burn=0,n_thin=1,
                     vcov_augment = 1,prog_bar = T,jumpyness = 1, seed = NA)
{

  if (is.na(seed) == F)
  {
    set.seed(seed)
  }

  par_curr <- par_init

  if (class(vcov) == 'data.frame'){ vcov <- as.matrix(vcov) }

  stream_pars = which(is.na(suppressWarnings(as.numeric(gsub('.*\\.', '',parnames)))) == F)

  par_mat <- matrix(NA,nrow = 20,ncol = 3)

  par_mat[,1] <- 1:20

  par_mat[,2:3] <- (par_init[,stream_pars])

  colnames(par_mat) <- c('stream','mu','alpha')

  ll_curr <- mozzy_like(pars = par_curr, par_mat = par_mat, dat = dat,parnames = parnames, stream_pars = stream_pars)$loglike

  outvars <- c(tolower(colnames(par_curr)),'ll')

  par_out <- matrix(NA, nrow = (n_sim - n_burn), ncol = length(outvars))

  i_count <- 0

  if(prog_bar == T)
  {
    pb <- txtProgressBar(min = 1, max = n_sim, style = 3)
  }

  #   a <- proc.time()

  orig_vcov <- vcov

  sel_count <- 0

  for (i in 1 : n_sim)
  {

    if (i < n_burn) {
      vcov <- orig_vcov * vcov_augment
    }else{
      vcov <- orig_vcov
    }

    par_next <- rmvnorm(n = 1, mean= par_curr, sigma =  vcov)

    colnames(par_next) <- parnames

#     if(par_next[,'mu.f'] < 0 ){browser()}

    par_mat[,2:3] <- (par_next[,stream_pars])

    ll_next <- mozzy_like(pars = par_next, par_mat = par_mat,dat = dat,parnames = parnames,stream_pars = stream_pars)$loglike

        if(prog_bar == T)
        {

    setTxtProgressBar(pb, i)

        }

    rand_accept <- log(runif(1,0,jumpyness))

    if (ll_next > ll_curr + rand_accept)
    {
      par_curr <- par_next
      ll_curr <- ll_next
      sel_count <- sel_count + 1
    }

    #     if (i > n_burn & i %% n_thin == 0) #%% means divisible by
    if (i > n_burn) #%% means divisible by
    {
      i_count <- i_count + 1
      par_out[i_count,] <- c(as.matrix(par_curr),ll_curr)
    }
  }
  # proc.time() - a
  par_out <- as.data.frame(par_out)

  colnames(par_out) <- outvars

  accepted_runs <- data.frame(num_selected = sel_count, num_simulated = n_sim, perc_selected = sel_count/n_sim,
                              perc_rejected = 1 - sel_count/n_sim, num_kept = n_sim - n_burn)

  return(list(posteriors = par_out[1:i_count,], initial_par = par_init,
              accepted_runs = accepted_runs ))
}