#' A function to return the likelihood for assignmen
#'
#' \code{mosquito.like} returns the negative log likelood of the data
#' @param pars is a list of parameters for the model
#' @param dat the stream-distance-count data
#' @param parnames vector of parameter names
#' @param stream_pars index of parameters that are stream specific
#' @param dat data
#' @param sigma.hyperprior.mu

mozzy_like <- function(pars ,parnames, par_mat,stream_pars, dat , sigma.hyperprior.mu = 1,sigma.hyperprior.sd = 0.1,
                       alpha.hyperprior.mu = .01, alpha.hyperprior.sd = 0.1)
{

  # Data likelihood ----

  pred_dat <- dat

  pred_dat$pred_count <- par_mat[dat$stream,'mu'] * exp(par_mat[dat$stream,'alpha'] * dat$distance)

  pred_dat$ll <- dpois(pred_dat$count,pmax(0,pred_dat$pred_count), log = T)

  #   pred_plot <- pred_dat %>%
  #   ggplot(aes(distance,count)) +
  #     geom_point() +
  #     geom_line(aes(distance, pred_count)) +
  #     facet_wrap(~stream)

  # Mu prior ----

  mu_expected <- pars[,'mu.f'] + ( (pars[,'mu.l'] - pars[,'mu.f']) * (par_mat[,'stream']-1) ) / (20-1)

  mu_ll <- dnorm(par_mat[,'mu'], mean = mu_expected, sd = pars[,'sigma.mu'], log = T)

  # Alpha prior ----

  alpha_ll <- dnorm(par_mat[,'alpha'],mean = pars[,'alpha.bar'], sd = pars[,'sigma.alpha'], log = T)

  # Hyperprior likelihoods ----

  mu_f_ll <- dunif(pars[,'mu.f'],0,1e10, log = T)

  mu_l_ll <- dunif(pars[,'mu.l'],0,1e10, log = T)

  sigma_mu_ll <- dnorm(pars[,'sigma.mu'], sigma.hyperprior.mu,sigma.hyperprior.sd, log = T)

  sigma_alpha_ll <- dnorm(pars[,'sigma.alpha'], alpha.hyperprior.mu,alpha.hyperprior.sd, log = T)

  ll <- sum(pred_dat$ll) + sum(mu_ll) + sum(alpha_ll) +
    sigma_mu_ll + sigma_alpha_ll   + mu_f_ll + mu_l_ll

  return(list(loglike = ll,data_loglike = sum(pred_dat$ll), deviance = -2*sum(pred_dat$ll), pred_dat = pred_dat))

}