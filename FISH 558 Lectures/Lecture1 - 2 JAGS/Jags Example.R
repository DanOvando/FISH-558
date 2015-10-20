
# JAGS example ------------------------------------------------------------



# Read in Data ------------------------------------------------------------

rm(list = ls())
library(runjags)
library(coda)
library(ggmcmc)
ogdat <- read.table('Lect2.dat', header = T)

head(ogdat)

summary(ogdat)

n <- dim(ogdat)[1]

dat <- list(plot = ogdat$Plot, effort = ogdat$Effort, count = ogdat$Count,n = n,
            numpars = 10)

inits <- list(lambda=rep(log(5),10), mean = rep(5,10), tau = rep(10,10),
              .RNG.name="base::Super-Duper", .RNG.seed=1, log.hypermu = 5, precision = 1)

model <- 'model{
for (i in 1:n)
{
  mu[i] <- lambda[plot[i]] * effort[i] #Distribution for estimating the predicted counts
  count[i] ~ dpois(mu[i]) #predicted counts
}
for (p in 1:numpars)
{
lambda[p] ~ dlnorm(log.hypermu,precision) #by patch catcheability thing
}
log.hypermu ~ dunif(-1000,1000)
precision ~ dunif(0,1000.0)
hypermu <- exp(log.hypermu)
sigma <- 1/sqrt(precision)
  }'


results <- run.jags(model=model, monitor=c("lambda",'log.hypermu','precision','hypermu','sigma','mu'),
                    data=dat, n.chains=1, method="rjags", inits=inits,plots=T,monitor.deviance=T,
                    silent.jag=F,burnin = 10000,thin = 10)

# results <- run.jags(model=model, monitor=c('mu'),
#                     data=dat, n.chains=1, method="rjags", inits=inits,plots=T,monitor.deviance=T,
#                     silent.jag=F,burnin = 10000,thin = 10)


# Predict out -------------------------------------------------------------

mcmc <- data.frame(as.matrix(results$mcmc))

head(mcmc)


lambda = gather(as.data.frame(((mcmc[,grepl('lambda.',colnames(mcmc))])),stringsAsFactors = F),convert = T)

real <- data.frame(plot = dat$plot,effort = dat$effort,real = dat$count,key = paste('lambda.',dat$plot,'.', sep = ''),stringsAsFactors = F)

real$data.point <- 1:dim(real)[1]

count <- left_join(lambda,real,by = 'key')

rands <- rpois(n = n * dim(mcmc)[1], lambda = count$value * count$effort ) #predict catches

posterior_preditive_counts <- data.frame(count,predicted.count = rands)
#
# mu_posterior_predictive <- gather(data.frame(
#   matrix(nrow = dim(mcmc)[1],ncol = n,rands,byrow = F)),convert = T)
#

# groups <- post.predict %>%
#   group_by(key) %>%
#   summarize(mv = max(value), mr = mean(real))

# quartz()
posterior.predictive.plot <- ggplot(subset(posterior_preditive_counts), aes(x = predicted.count,fill = data.point)) +
  geom_histogram() +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~data.point, scales = 'free_x')+ theme(legend.position = 'none', axis.text.x= element_text(size = 6))
ggsave('posterior_predictive_plot.pdf', plot = posterior.predictive.plot)

# quartz()
posterior.predictive.vplot <- (ggplot(posterior_preditive_counts, aes(factor(data.point),predicted.count)) + geom_violin() +
                                 geom_point(data = real,aes(x = factor(data.point), y = (real)), size = 2, shape = 21, fill = 'lightseagreen') +
                                 theme(axis.text.x = element_text(size = 5)) + xlab('Data Point') + ylab('Counts'))
ggsave('posterior_predictive_vplot.pdf', plot = posterior.predictive.vplot)



# ggmcmc(ggs(results$mcmc),file="test jags fit.pdf")
#
#   print(results)
# #   plot(results,type="trace",layout=c(2,2))
#
#   # quartz()
#   plot(results, layout = runjags.getOption("plot.layout"),
#        new.windows = runjags.getOption("new.windows"), file = "testjags.pdf")

#   crosscorr.plot(results)







