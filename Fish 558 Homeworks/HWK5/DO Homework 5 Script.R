
set.seed(54321)

run <- '1.2'

num_cores <- 4 #Number of cores to run in parallel

run_folder <- paste('Results/',run,'/',sep = '')

run_sims <- F # Set to false if you want to run off a stored run

run_on_server <- F

if (run_on_server == T){run_sims <- T}

if (dir.exists(run_folder) == F){dir.create(run_folder, recursive = T)}

library(knitr)
knitr::opts_chunk$set(fig.path='Figs/', echo=FALSE, warning=FALSE, message=FALSE)
library(gridExtra)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(broom)
library(foreach)
library(scales)
library(stargazer)
library(mvtnorm)
library(proftools)
library(shiny)
library(doMC)
library(DT)


#' Convert moment to shape parameters
#'
#' \code{moments_to_shape} takes the mean and cv
#' of a beta distribution and converts it to the shape
#' parameters
#' @param mu the mean of the beta distribution
#' @param cv the cv of the beta distribution

moments_to_shape <- function(mu,cv)
{

  var <- (mu*cv)^2 #calculate variance

  a <- (((1/mu - 1)*mu^3)/var) - mu

  b <- a*(1/mu - 1)

  #   a <- ( (1-mu)/var - 1/mu )*mu^2
  #
  #   b <- a*(1/mu - 1)

  if (mu >=1){warning('mu is greater than or equal to 1, reduce it!')}

  return(list(shape1 = a, shape2 = b))
}

#' runs problembeast model
#'
#'\code{sim_problembeast} simulates populations of
#'the endangered problembeast forward under
#'different states of nature and hunting
#'@param pb a list of problembeast population
#'parameters
#'@param sim_years number of years to run
#'@param hunted_patch the patch in which hunting
#'occurs
#'@param s_j juvenile survival
#'@param _value the price per hunted animal
#'at various ages
#'@param results blank matrix to store results
#'@param recruits index of recruits
#'@param juveniles index of juveniles
#'@param adults index of all adults
#'@param mature index of mature adults
#'@param plus index of the plus group
#'
sim_problembeast <- function(pb,sim_years = 101,hunted_patch = 0,s_j = .75,
                             recruit_value = 1, juv_value = 5, adult_value = 30,
                             results, recruits,juveniles,adults,mature,plus, animals)
{

  pb$s_j <- s_j #input juvenile survival

  results$year <- 1:(sim_years+1) #input yeats

  results$hunted_patch <- hunted_patch #mark hunted patch

  results$s_j <- s_j #store juvenile survival

  phi <- matrix(1, nrow = 3, ncol = sim_years)

  if (hunted_patch>0){ # If there is any hunting generate a vector of hunting survivals over time for the appropriate patch

    shape1 <- eval(parse(text = paste('pb$phi_',hunted_patch,'_shape',sep = '') ))$shape1

    shape2 <- eval(parse(text = paste('pb$phi_',hunted_patch,'_shape',sep = '') ))$shape2

    phi[hunted_patch,] <-  rbeta(sim_years,shape1 = shape1, shape2 = shape2)


  }

  #Store initial conditions

  n_adults <- sum(pb$n_at_age[1, mature])

  results$adults[1] <- n_adults

  results$animals[1] <- sum(pb$n_at_age[1, animals])

  pb$n_at_age <- as.matrix(pb$n_at_age)

  for (y in 1:(sim_years)) #loop away
  {

    # Juvenile component ----

    juv_survive <- rbinom(length(juveniles),(pb$n_at_age[y,juveniles - 1]),pb$s_j)

    juv_survive_hunting <- rbinom(length(juveniles),(juv_survive),phi[2,y])

    juv_hunted <- juv_survive - juv_survive_hunting

    pb$n_at_age[y + 1,juveniles] <- juv_survive_hunting

    # Adult component ----

    adults_survive <- rbinom(length(adults),(pb$n_at_age[y,adults - 1]),pb$s_a)

    adults_survive_hunting <- rbinom(length(adults), (adults_survive),phi[3,y])

    adults_hunted <- adults_survive - adults_survive_hunting

    pb$n_at_age[y + 1, adults] <- adults_survive_hunting

    # Plus group component ---

    enter_plus_survive <- rbinom(length(plus),(pb$n_at_age[y,plus - 1]),pb$s_a)

    enter_plus_survive_hunting <- rbinom(length(plus), (enter_plus_survive),phi[3,y])

    enter_plus_hunted <- enter_plus_survive - enter_plus_survive_hunting

    enter_plus <- enter_plus_survive_hunting

    stay_plus_survive <- rbinom(length(plus),(pb$n_at_age[y,plus]),pb$s_a)

    stay_plus_survive_hunted <- rbinom(length(plus), (stay_plus_survive),phi[3,y])

    stay_plus_hunted <- stay_plus_survive - stay_plus_survive_hunted

    stay_plus <- stay_plus_survive_hunted

    pb$n_at_age[y + 1, plus] <- enter_plus + stay_plus

    n_adults <- sum(pb$n_at_age[y + 1, mature])

    # Recruitment comonent ----

    recruits_born <- rbinom(length(recruits), n_adults, pb$b)

    recruits_survive_hunting <- rbinom(length(recruits), recruits_born,phi[1,y])

    recruits_hunted <- recruits_born - recruits_survive_hunting

    pb$n_at_age[y + 1,recruits] <- recruits_survive_hunting

    results$adults[y+1] <- n_adults

    # Process hunting ----

    hunted_in_two <- sum(juv_hunted)

    hunted_in_three <- sum(adults_hunted) + sum(enter_plus_hunted) + sum(stay_plus_hunted)

    hunted <-  hunted_in_two + hunted_in_three

    hunted_value <- hunted_in_two * juv_value + hunted_in_three * adult_value

    if (hunted_patch >1){
      results$hunted[y] <- hunted
      results$hunt_value[y] <- hunted_value
    }
    if(hunted_patch == 1)
    {
      results$hunted[y + 1] <- recruits_hunted
      results$hunt_value[y+1] <- recruits_hunted * recruit_value
    }

    results$animals[y+1] <- sum(pb$n_at_age[y + 1,animals])

    results$plus_group[y+1] <- pb$n_at_age[y + 1,plus]
  } #close year loop

  results <- results[1:sim_years,]

  return(results)

}

sim_years <- 101

ages <- 16

n_sims <- 1000

problembeast <- list()

problembeast$n_at_age <- as.data.frame(matrix(NA,nrow = sim_years+1, ncol = ages + 1))

colnames(problembeast$n_at_age) <- c('year',paste('age.',0:(ages-1), sep = ''))

problembeast$n_at_age[1,] <- c(1,98, 68, 48, 38, 27, 24, 16, 11,
                               12, 13, 3, 11, 1, 5, 2, 113)

problembeast$s_j <- NA

problembeast$s_a <- 0.95

problembeast$b <- 0.75

problembeast$phi_1_mean <- 0.7

problembeast$phi_1_cv <- 0.1

problembeast$phi_1_shape <- moments_to_shape(mu = problembeast$phi_1_mean,cv = problembeast$phi_1_cv)


problembeast$phi_2_mean <- 0.95

problembeast$phi_2_cv <- 0.1

problembeast$phi_2_shape <- moments_to_shape(mu = problembeast$phi_2_mean,cv = problembeast$phi_2_cv)

problembeast$phi_3_mean <- 0.97

problembeast$phi_3_cv <- 0.1

problembeast$phi_3_shape <- moments_to_shape(mu = problembeast$phi_3_mean,cv = problembeast$phi_3_cv)


resultnames <- c('year','adults','animals','hunted','hunt_value','hunted_patch','s_j','plus_group')

resultmat <- as.data.frame(matrix(0,nrow = sim_years+1,ncol = length(resultnames)))

colnames(resultmat) <- resultnames

recruits <- which(colnames(problembeast$n_at_age) %in% 'age.0')

animals <- (grep('age.',colnames(problembeast$n_at_age)))

juveniles <- which(colnames(problembeast$n_at_age) %in% paste('age',1:8,sep = '.'))

adults <- which(colnames(problembeast$n_at_age) %in% paste('age',9:14,sep = '.'))

mature <- which(colnames(problembeast$n_at_age) %in% paste('age',9:15,sep = '.'))

plus <- which(colnames(problembeast$n_at_age) == 'age.15')


sims <- expand.grid(s_j = c(0.75,0.8,0.82),hunted_patch = c(0,1,2,3),sim = 1:n_sims)

run_sim <- function(j,sims,problembeast,resultmat,recruits,juveniles,adults,plus,mature)
{
  out <- sim_problembeast(pb = problembeast,
                          hunted_patch = sims$hunted_patch[j],s_j = sims$s_j[j], results = resultmat,juveniles = juveniles,adults = adults, mature = mature, plus = plus,recruits = recruits,
                          animals = animals)

  out$sim <- j


  return(out)
}

if (run_sims == T){
  registerDoMC(cores = num_cores)


  results <- (foreach(i = 1:dim(sims)[1]) %dopar%
                (run_sim(i,sims = sims,problembeast = problembeast, resultmat = resultmat,juveniles = juveniles,
                         adults = adults, mature = mature, plus = plus, recruits = recruits)))

  if (run_on_server == F)
  {
    save(file = paste(run_folder,'problembeast simulations.Rdata',sep = ''), results)
  }
}else{
  if (run_on_server == F)
  {
    load(file = paste(run_folder,'problembeast simulations.Rdata',sep = ''))
  }
}

tidy_results <- ldply(results) %>%
  group_by(sim) %>%
  mutate(cumu_hunt_value = cumsum(hunt_value)) %>%
  ungroup()



decision_table <- tidy_results %>%
  group_by(sim) %>%
  mutate(cumu_hunt_value = cumsum(hunt_value)) %>%
  ungroup() %>%
  subset(year == 101) %>%
  group_by(s_j,hunted_patch) %>%
  summarize(mean_adults = mean(adults),p_lessthan_1000 = mean(adults < 1000), mean_value = mean(cumu_hunt_value), mean_animals = mean(animals)) %>%
  ungroup() %>%
  group_by(hunted_patch) %>%
  mutate(expected_adults = mean(mean_adults), expected_p_lt_1000 = mean(p_lessthan_1000),
         expected_value = mean(mean_value), expected_animals = mean(mean_animals)) %>%
  ungroup() %>%
  mutate(state_of_nature = paste('juv survival = ',s_j, sep = ''))

decision_table1 <- select(decision_table, state_of_nature,hunted_patch,mean_adults,expected_adults) %>%
  spread(state_of_nature,mean_adults) %>%
  mutate('Expected Adults' = expected_adults) %>%
  select(-expected_adults) %>%
  rename('Hunted Patch' = hunted_patch)

decision_table2 <- select(decision_table, state_of_nature,hunted_patch,p_lessthan_1000,expected_p_lt_1000) %>%
  spread(state_of_nature,p_lessthan_1000) %>%
  mutate('Expected p(A < 1000)' = expected_p_lt_1000) %>%
  select(-expected_p_lt_1000) %>%
  rename('Hunted Patch' = hunted_patch)

decision_table3 <- select(decision_table, state_of_nature,hunted_patch,mean_animals,expected_animals) %>%
  spread(state_of_nature,mean_animals) %>%
  mutate('Expected Animals' = expected_animals) %>%
  select(-expected_animals) %>%
  rename('Hunted Patch' = hunted_patch)

decision_table_value <- select(decision_table, state_of_nature,hunted_patch,mean_value,expected_value) %>%
  spread(state_of_nature,mean_value) %>%
  mutate('Expected Hunting Value' = expected_value) %>%
  select(-expected_value) %>%
  rename('Hunted Patch' = hunted_patch)

tradeoff <- (ggplot(decision_table,aes(mean_adults,mean_value, shape = factor(hunted_patch),color = s_j)) +
               geom_point(size = 4, alpha = 0.75) +
               #                geom_line(aes(mean_adults,mean_value,color = s_j)) +
               scale_color_gradient(low = 'red',high = 'green') +
               theme_light())

tradeoff2 <- (ggplot(subset(decision_table, hunted_patch >0),aes(1-p_lessthan_1000,mean_value,fill = factor(hunted_patch))) +
                geom_point(shape = 21, size = 4, alpha = 0.75) +
                scale_fill_brewer(name = 'Hunted Patch',palette = 'Spectral') +
                facet_grid(state_of_nature~., scales = 'free') +
                theme_light() +
                xlab('P(Adults > 1000)') +
                ylab('Value of Hunting'))

tradeoff3 <- (ggplot(subset(decision_table, hunted_patch >0),aes(1-expected_p_lt_1000,expected_value,fill = factor(hunted_patch))) +
                geom_point(shape = 21, size = 4, alpha = 0.75) +
                scale_fill_brewer(name = 'Hunted Patch',palette = 'Spectral') +
                #                 facet_grid(state_of_nature~., scales = 'free') +
                theme_light() +
                xlab('P(Adults > 1000)') +
                ylab('Value of Hunting'))

tidy_results$hunted_patch_name <- paste('Hunted Patch = ',tidy_results$hunted_patch, sep = '')

outcome_plot <- subset(tidy_results, year == 101) %>%
  ggplot(aes(factor(hunted_patch),adults, fill = factor(s_j))) +
  geom_boxplot(alpha = 0.75) +
  xlab('Hunted Patch') +
  ylab('Number of Adults') +
  scale_fill_brewer(name = 'Juv. Survival',palette = 'Spectral') +
  theme_light() +
  theme(legend.position = 'top') +
  geom_hline(aes(yintercept = 1000), linetype = 'dashed')

outcome_plot_value <- subset(tidy_results, year == 101) %>%
  subset(hunted_patch >0) %>%
  ggplot(aes(factor(hunted_patch),cumu_hunt_value, fill = factor(s_j))) +
  geom_boxplot(alpha = 0.75) +
  xlab('Hunted Patch') +
  ylab('Cumulative Hunting Value') +
  scale_fill_brewer(name = 'Juv. Survival',palette = 'Spectral') +
  theme_light() +
  theme(legend.position = 'top')

with_information <- decision_table %>%
  subset(hunted_patch > 0) %>%
  group_by(s_j) %>%
  summarise(best_strat = hunted_patch[mean_animals == max(mean_animals)],
            val_at_best = mean_value[mean_animals == max(mean_animals)],
            best_animals = mean_animals[mean_animals == max(mean_animals)])

expected_animals_with_info <- mean(with_information$best_animals)

expected_value_with_info <- mean(with_information$val_at_best)

without_information <- decision_table %>%
  subset(hunted_patch >0) %>%
  subset(expected_animals == max(expected_animals))

expected_animals_wout_info <- mean(without_information$mean_animals)

expected_value_wout_info <- mean(without_information$mean_value)

voi <- expected_animals_with_info - expected_animals_wout_info
