# Wrapper script to run in terminal ----
rm(list = ls())

run <- '3.3'

sims <- 1.5e6

run_folder <- paste('Results/',run,'/',sep = '')

if (dir.exists(run_folder) == F){dir.create(run_folder, recursive = T)}


# set.seed(54321)
library(knitr)
# knitr::opts_chunk$set(fig.path='Figs/', echo=FALSE, warning=FALSE, message=FALSE)
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

multichain_mcmc_results <- multi_mozzy_mcmc(par_guess = as.matrix(par_guess),jitfactor = 5,
                                            num_starts = 3,numcores = 3, parnames = colnames(par_guess),
                                            dat = dat, vcov = vcov, prog_bar = F, n_sim = sims,
                                            n_burn = round(.6*sims,0), n_thin = 1,
                                            vcov_augment = (2.4/sqrt(45))^2,jumpyness = .001)

print('finished parallel chains')

save(file = paste(run_folder,'multi_chain_mcmc.Rdata',sep = ''), multichain_mcmc_results)
