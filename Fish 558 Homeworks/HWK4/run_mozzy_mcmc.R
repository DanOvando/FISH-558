# Wrapper script to run in Terminal

#Initials ----

sims <- 1.5e6

run <- '3.31'

run_folder <- paste('Results/',run,'/',sep = '')

if (dir.exists(run_folder) == F){dir.create(run_folder, recursive = T)}

set.seed(54321)

library(knitr)
knitr::opts_chunk$set(fig.path='Figs/', echo=FALSE, warning=FALSE, message=FALSE)
library(gridExtra)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(broom)
library(coda)
library(ggmcmc)
# library(LaplacesDemon)
library(foreach)
library(scales)
library(stargazer)
library(mvtnorm)
library(doMC)
library(proftools)
source('mozzy_likelihood.R')
source('mozzy_mcmc.R')
source('multi_mozzy_mcmc.R')
source('thin_mcmc.R')

# Data ----

dat <- read.csv(file = 'hwk4_data.csv',stringsAsFactors = F) %>%
  dplyr::rename(stream = Steam) %>%
  gather('distance','count', grep('X',colnames(.))) %>%
  mutate(distance = as.numeric(gsub('X','',distance)))

colnames(dat) <- tolower(colnames(dat))

par_guess <- (read.csv(file = 'hwk4_pars.csv', stringsAsFactors = F))

colnames(par_guess) <- tolower(colnames(par_guess))


vcov <- read.csv(file = 'hwk4_vcov.csv', stringsAsFactors = F)

rownames(vcov) <- colnames(vcov)


dat_plot <- dat %>%
  ggplot(aes(distance, count, fill = stream, size = count)) +
  geom_point(shape = 21) +
  facet_wrap(~stream)

vcov.plot <- vcov %>%
  mutate(var1 = colnames(.)) %>%
  gather('var2','covar',which(grepl('var1',colnames(.)) == F)) %>%
  ggplot(aes(var1,var2,fill = covar)) +
  geom_tile()

# MCMC ----


#
# optim(.3,tune_scalar,n_sim = n_sim,par_init = as.matrix(par_guess),parnames = colnames(par_guess),dat = dat,vcov = vcov, n_burn=round(.5*n_sim,0)
#       ,n_thin=1,prog_bar = F,jumpyness = 1, seed = NA,targ_accept_rate = 0.3, lower = .05, upper = 1.5)
#
# tune_scalar(.3,n_sim = 1000,par_init = as.matrix(par_guess),parnames = colnames(par_guess),dat = dat,vcov = vcov, n_burn=round(.5*n_sim,0)
# ,n_thin=1,prog_bar = F,jumpyness = 1, seed = NA,targ_accept_rate = 0.3)

# a <- proc.time()

mcmc_results <- mozzy_mcmc(par_init = as.matrix(par_guess),
                           parnames = colnames(par_guess), dat = dat,
                           vcov = vcov, prog_bar = T, n_sim = sims, n_burn = round(.6*sims,0),
                           n_thin = 1,vcov_augment = (2.4/sqrt(45))^2,jumpyness = 1,targ_accept_rate = .2)

# proc.time() - a


mcmc_posteriors <- thin_mcmc(chains = mcmc_results$posteriors, thin_every = sims/1000)

ggmcmc(ggs(mcmc(mcmc_posteriors)), file=paste(run_folder,"mozzy_mcmc_diagnostics.pdf", sep = ''))


save(file = paste(run_folder,'mozzy_mcmc.Rdata', sep = ''), mcmc_results)


