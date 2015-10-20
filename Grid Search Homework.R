#Homework Grid Search Methods
#

rm(list = ls())
set.seed(100)
library(gridExtra,quietly = T)
library(ggplot2, quietly = T)
library(plyr, quietly = T)
library(dplyr, quietly = T)
library(tidyr, quietly = T)
library(knitr, quietly = T)

like_fun <- function(pred,obs,sigma,n)
{
#   nll <- n * log(sigma) + sum( (log(pred) - log(obs))^2 / 2*sigma^2)
  nll <- sum((log(pred) - log(obs))^2)
}


dat <- read.table(file = 'Fish 558 Workshop/Workshop Day 3/GRAY.TXT',header = F,stringsAsFactors = F)

colnames(dat) <- c('year','pop')

summary(dat)

ggplot(dat,aes(year,pop)) + geom_point() + xlab('Year') + ylab('Whale Population')

sigma <- 0.1

# Part 1

grid_search <- function(grid_space = 5,lower_pgrid = 10000,upper_pgrid = 15000,lower_alpha = 0,upper_alpha = 0.05)
{

pgrid = seq(lower_pgrid,upper_pgrid,length.out = grid_space)

alpha_grid <- seq(lower_alpha,upper_alpha,length.out = grid_space)

years <- dat$year

grid <- expand.grid(p1968 = pgrid,alpha = alpha_grid,year = years)

grid$pophat <- with(grid,p1968 * exp(alpha * (year - 1968)) * exp(rnorm(1,0,.1)))

grid <- join(grid,dat,by = 'year')

grid$likelihood <- with(grid,exp(-(log(pophat) - log(pop) )^2 ))
browser()
grid <- grid %>%
  group_by(alpha,p1968) %>%
  summarize(total_like = prod(likelihood)) %>%
  ungroup() %>%
  mutate(prior = dnorm(alpha,mean = .04,sd = .01) * dnorm(p1968,mean = 12000,sd = 1000)) %>%
  mutate(posterior = prior * total_like, norm_posterior = posterior/sum(posterior))


alpha_like <- grid %>%
  group_by(alpha) %>%
  summarize(marg_like = sum(norm_posterior))

p1968_like <- grid %>%
  group_by(p1968) %>%
  summarize(marg_like = sum(norm_posterior))

return(list(grid = grid, alpha_like = alpha_like,p1968_like = p1968_like))
}

grid_space <- 100

grid <- grid_search(grid_space = grid_space)

quartz()
ggplot(grid$alpha_like,aes(alpha,marg_like)) + geom_point()

quartz()

ggplot(grid$p1968_like,aes(p1968,marg_like)) + geom_point()

quartz()

persp(x = unique(grid$grid$alpha), y = unique(grid$grid$p1968), z = matrix(grid$grid$norm_posterior,nrow = grid_space, ncol = grid_space))

# Part 2







