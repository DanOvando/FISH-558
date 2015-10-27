## Wrapper script to reproduce Bowhead Model outside of markdown

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
library(LaplacesDemon)
library(foreach)
library(scales)
devtools::load_all('BowheadModel')

set.seed(12345)

runname <- 'Version 1.0'

runfolder <- paste('Results',runname, sep = '/')

if (dir.exists(runfolder) == F) {dir.create(runfolder, recursive = T)}

# Load in and Process Data ------------------------------------------------------------

catch.dat <- read.csv('Hwk3.csv',stringsAsFactors = F, header = F)
colnames(catch.dat) <- c('year','catch')

abund.dat <- read.csv('HWK3B.csv',stringsAsFactors = F, header = F)
colnames(abund.dat) <- c('year','abundance', 'abundance.cv')

dat <- full_join(catch.dat,abund.dat, by = 'year') %>%
  subset(is.na(year) == F)

whale.data.plot <- dat %>%
  gather('metric','whale.numbers',2:3) %>%
  ggplot(aes(year,whale.numbers, fill = metric)) +
  geom_point(shape = 21)

a <- proc.time()

# Make object with all the stuff for this particular whale population simulation
whales <- make.whales(dat = dat,catch.dat = catch.dat,f.max = .34,s.0 = .2,s.rest = .97, extra.time = NA,K = 10000, use.catch = F)

#Run the model using the created object
whale.proj <- whale.pop.model(whales = whales, dat = dat, n.years = dim(dat)[1])

# Calculate the negative log likelihood
nll <- whale.likelihood(whale.proj)

# nlminb(start = c(.9,.9,2000,.3),whale.optim, dat = dat, catch.dat= catch.dat, lower = c(.8,.8,1000,.1), upper = c(1,1,3000,.4))

proc.time() - a

#Check whether population is stable at carrying capacity
b <- whale.proj$pop %>%
  dplyr::select(-numbers) %>%
  gather('age','numbers',age.0:age.13) %>%
  subset(age != 'age.0') %>%
  group_by(year) %>%
  summarize(total.whales = sum(numbers), mature.whales = last(numbers)) %>%
  ggplot(aes(year,total.whales)) +
  geom_point() +
  #   ylim(0,NA) +
  geom_hline(aes(yintercept = (whales$life$K)))



# Run post-model pre-data SIR ---------------------------------------------

binary.whales <- whale.sir(Nout = 1000, k.lower = 10000, dat = dat, catch.dat = catch.dat, progbar = T, mode = 'simple', MLE = log(0.99), use.catch = T)

# Run post-data SIR ---------------------------------------------

functional.whales <- whale.sir(Nout = 200, k.lower = 10000, dat = dat, catch.dat = catch.dat, progbar = T, mode = 'fit', MLE = 2.22)

# Simulate whale populations ----------------------------------------------

future.catches <- c(67, 134, 201)

future.whales <- foreach(f = 1:length(future.catches), .combine = rbind) %do%
{
  simmed.whales <- sim.whales(dat = dat, catch.dat = catch.dat, possible.whales = functional.whales, extra.time = 21, extra.catch = future.catches[f])
}


# Process Results ---------------------------------------------------------

processed.binary.whales <- process.whale.SIR(fitted.whales = binary.whales, dat = dat, catch.dat = catch.dat, rungroup = 'Binary Whales', runfolder = runfolder, use.catch = T)


processed.functional.whales <- process.whale.SIR(fitted.whales = functional.whales, dat = dat, catch.dat = catch.dat, rungroup = 'Fitted Whales', runfolder = runfolder)

check.catches <- future.whales %>%
  group_by(year, catchlevel) %>%
  summarize(total.catch = mean(catch, na.rm = T))

catch.plot <- check.catches %>%
  ggplot(aes(year, total.catch, fill = factor(catchlevel))) +
  geom_point(shape = 21)


summarized.future.whales <- future.whales %>%
  group_by(year,catchlevel) %>%
  summarize(lower.5 = quantile(predicted.whales,probs = 0.05, na.rm = T),
            upper.95 =quantile(predicted.whales,probs = 0.95, na.rm = T),
            median.whales = median(predicted.whales)) %>%
  left_join(dplyr::select(dat, year,abundance), by = 'year')

summarized.future.whales$catchlevel[summarized.future.whales$year <= 2002] <- 'Historic'

future.whale.summary.plot <- summarized.future.whales %>%
  subset(year > 1975) %>%
  ggplot(aes(x = year)) +
  geom_ribbon(aes(ymin = lower.5, ymax = upper.95, color = factor(catchlevel),fill = factor(catchlevel)), alpha = 0.25) +
  geom_line(aes(year, median.whales, color = factor(catchlevel)), size = 1.5, linetype = 'longdash') +
  geom_point(aes(year, abundance), shape = 19, color = 'black') +
  scale_fill_manual(name = 'Catch Scenario', values =  c('lightgoldenrod','orangered','green4','steelblue2')) +
  scale_color_manual(name = 'Catch Scenario', values =  c('lightgoldenrod','orangered','green4','steelblue2')) +
  xlab('Year') +
  ylab('Numbers of Whales') +
  geom_vline(aes(xintercept = 2003), linetype = 'longdash')  +
  theme_bw()

ggsave(paste(runfolder,'Whale Projections.pdf',sep = '/'), plot = future.whale.summary.plot, height = 6, width = 8 )


whale.groups <- future.whales %>%
  subset(year == 2003) %>%
  mutate(state.of.nature = cut(predicted.whales, breaks = c(0,7000,8000,20000), labels = c('<7000','7000-8000','>8000')))

probs.of.nature <- whale.groups %>%
  group_by(state.of.nature) %>%
  summarize(n.occured = length(iteration)) %>%
  ungroup() %>%
  mutate(prob.occuring = n.occured / sum(n.occured))


# Construct decision table ------------------------------------------------

decision.table <- future.whales %>%
  left_join(dplyr::select(whale.groups,iteration,state.of.nature), by = 'iteration') %>%
  group_by(state.of.nature, catchlevel) %>%
  mutate(nvk =predicted.whales/K ) %>%
  summarize(prob.increased.whales = mean( (predicted.mature.whales[year == 2023] / predicted.mature.whales[year == 2003])>1),
            prob.nvk = mean((nvk)[year == max(year)] > 0.5 & nvk[2023]>nvk[2003])) %>%
  left_join(probs.of.nature, by = 'state.of.nature') %>%
  ungroup() %>%
  mutate(expected.increase = prob.increased.whales * prob.occuring,
         expected.nvk = prob.nvk * prob.occuring) %>%
  group_by(catchlevel) %>%
  mutate(expected.outcome = sum(expected.increase),
         expected.nvk.outcome = sum(expected.nvk)) %>%
  arrange(desc(expected.outcome))

save.image(file=paste(runfolder,'bowhead whale decision analysis.rdata',sep='/'))

write.csv(file=paste(runfolder,'bowhead whale decision table.csv',sep='/'), decision.table)

# Make some other plots ---------------------------------------------------

predat <- processed.binary.whales$posterior.distributions

predat$model <- 'Pre-Data'

postdat <- processed.functional.whales$posterior.distributions

postdat$model <- 'Post-Data'

pre.post.dat <- rbind(predat,postdat)

posterior.distributions.plot  <- pre.post.dat %>%
  ggplot(aes(value, fill = model)) +
  scale_fill_brewer(name = 'SIR Model', palette = 'Spectral') +
  #     geom_histogram(aes(y = ..density..) ,alpha = 0.6) +
  geom_density(alpha = 0.6) +
  facet_wrap(~variable, scales = 'free') +
  theme_bw() +
  xlab('Value')+
  ylab('Density')

posterior.distributions.plot

NvK.plot <- future.whales %>%
  subset(year == max(year)) %>%
  ggplot(aes(predicted.whales/K, fill = factor(catchlevel))) +
  #   geom_histogram(alpha = 0.6) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(name = 'Future Catches',values = c('green','yellow','red')) +
  theme_bw() +
  xlab('% of Unexploited Population') +
  ylab('Density') +
  scale_x_continuous(labels = percent) +
  geom_vline(aes(xintercept = .50))

NvK.plot

dt2 <- decision.table %>%
  select(catchlevel,prob.nvk,expected.nvk.outcome,state.of.nature) %>%
  rename('State of Nature' = state.of.nature, 'Future Catch' = catchlevel,
         'P(N2023/K > 0.5)' = prob.nvk, 'Expected Value' = expected.nvk.outcome)

