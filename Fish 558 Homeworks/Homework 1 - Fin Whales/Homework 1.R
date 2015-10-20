library(dplyr)
library(bbmle)
library(ggplot2)
library(tidyr)

catch.dat <- read.csv('Fish 558 Homeworks/hwk1_catch_data.csv',stringsAsFactors = F) %>%
  dplyr::select(year,catch)

head(catch.dat)

survey.dat <- read.csv('Fish 558 Homeworks/hwk1_survey_data.csv',stringsAsFactors = F)


whale.dat <- full_join(catch.dat,survey.dat,by = 'year')

head(whale.dat)

whale.plot <- whale.dat %>%
  gather('data','value',2:3)

mf_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="data") {
    value[value=="catch"] <- "Catch"
    value[value=="n"]   <- "Abundance Estimate"
  }
  return(value)
}

# ggplot(whale.plot,aes(year,value, fill = data)) +
#   geom_point(shape = 21) +
#   facet_grid(data ~.,labeller=mf_labeller) +
#   theme(legend.position = 'none')

ggplot(whale.plot,aes(year,value, color = data,shape = data)) +
  geom_point() +
  scale_color_discrete(name=element_blank(),
                       breaks=c("catch", "n"),
                       labels=c("Catch", "Abundance Estimate")) +
  scale_shape_discrete(name=element_blank(),
                       breaks=c("catch", "n"),
                       labels=c("Catch", "Abundance Estimate"))

fit.whale.model <- function(r,k,tau,b0,z,dat,use = 0){

  dat$nhat <- (pella.tom.model(r = r, k = k, b0 = k,z = 2.39,catch = dat$catch, time = dim(dat)[1]))

  dat$sigma <- sqrt(dat$u ^2 + tau ^ 2)

  dat$nll <- -log(1/dat$sigma) + ( log(dat$nhat) - log(dat$n) )^2 / (2 * dat$sigma^2)

  show(sum(-log(1/dat$sigma) + ( log(dat$nhat) - log(dat$n) )^2 / (2 * dat$sigma^2), na.rm = T))

#   show(prod((1/dat$sigma) * exp(- (log(dat$nhat) - log(dat$n))^2 / (2 * dat$sigma^2)), na.rm = T))

  if (use == '0')
  {
    output <- (as.numeric(sum(dat$nll, na.rm = T)))
  }
  if (use != '0'){
    output <- dat
  }

  return(output)
}

fit.whale.model2 <- function(pars,b0 = NA,z = 2.39,dat,use = 'optim'){

  r <- pars[1]
  k <- pars[2]
  tau <- pars[3]

  dat$nhat <- (pella.tom.model(r = r, k = k, b0 = k,z = z,catch = dat$catch, time = dim(dat)[1]))

  dat$sigma <- sqrt(dat$u ^2 + tau ^ 2)

  dat$nll <- -log(1/dat$sigma) + ( log(dat$nhat) - log(dat$n) )^2 / (2 * dat$sigma^2)

  if (use == 'optim')
  {
    output <- -sum(dnorm(dat$n,mean = dat$nhat, sd = dat$sigma, log = T), na.rm = T)

#     output <- (as.numeric(sum(dat$nll, na.rm = T)))
  }
  if (use != 'optim'){
    output <- dat
  }

  return(output)
}


pella.tom.model <- function(r,k,b0,z,catch,time){

  b <-   matrix(NA,nrow = time,ncol = 1)

  b[1] <- b0

  for (t in 2:time)
  {
    last.b <- b[t-1]

    b[t] <- pmax(last.b + last.b * r * (1- (last.b/k)^z) - catch[t-1],0.0001)
  }

  return(b)
}

whale.fit2 <- nlminb(start = c(.1,1.2*max(whale.dat$n, na.rm = T),25),fit.whale.model2,dat = whale.dat)

whale.fit = mle(fit.whale.model,start = list(r = 0.1, k = ( 1.1*max(whale.dat$n, na.rm = T)), tau = (0.2)),
    fixed = list(z = 2.39,b0 = NA,dat = whale.dat), lower = c(.000001,1000,.0001),upper = c(2,20000,100))

fitted.whale <- fit.whale.model(r = coef(whale.fit)["r"],k = coef(whale.fit)["k"],
                                tau = coef(whale.fit)["tau"],dat = whale.dat,z = 2.39,use = 1,b0 = coef(whale.fit)["k"])

fitted.whale <- fit.whale.model2(pars = whale.fit2$par,dat = whale.dat,use = 'model')

ggplot(fitted.whale,aes(year,n)) +
  geom_point(aes(fill = 'Observed')) +
  geom_line(aes(year,nhat,color = 'Predicted')) +
  theme(legend.title = element_blank()) +
  xlab('Year') +
  ylab('Number of Fin Whales')



