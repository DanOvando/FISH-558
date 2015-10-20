schaefer_grid_search <- function(grid_space = 5,lower_r = .15,upper_r = .45,lower_k = 5000,upper_k = 8000,
                                 lower_phi = .8,upper_phi = 1.2,sigma = 0.1,dat)
{

  r_grid = seq(lower_r,upper_r,length.out = grid_space)

  k_grid <- seq(lower_k,upper_k,length.out = grid_space)

#   phi_grid <- seq(lower_phi,upper_phi,length.out = grid_space)

  years <- dat$year

  grid <- expand.grid(r = r_grid,k = k_grid,year = years) # create matrix of factorial combinations of alpha and o1968

  call.pop <- function(i,grid,dat,time)
  {
    pars <- log(c(grid$r[i], grid$k[i],1))
    outpop <- schaefer.model(pars,dat,time,mode = 'pop')
    return(outpop)
  }

  Pops <- lapply(1:dim(grid)[1],call.pop,grid = grid, dat = dat,time = dim(dat)[1]) %>% ldply()

  #   grid$pophat <- with(grid,p1968 * exp(alpha * (year - 1968))) #predict whale population in each year

  #   grid <- join(grid,dat,by = 'year') #add in observed data

  grid <- Pops %>%
    ungroup() %>%
    group_by(r,k) %>%
    summarize(total_like = mean(total_like)) %>% #Get total likelihood across all years for each alpha and p1968
    ungroup() %>%
    #     mutate(prior = dnorm(alpha,mean = .04,sd = .01) * dnorm(p1968,mean = 12000,sd = 1000)) %>% #calculate prior
    #     mutate(posterior = prior * total_like, norm_posterior = posterior/sum(posterior)) #calculate the normlized posterior
    mutate(posterior = total_like, norm_posterior = posterior/sum(posterior)) #calculate the normlized posterior
  # Marginal likelihood of alpha
  r_like <- grid %>%
    group_by(r) %>%
    summarize(marg_like = sum(norm_posterior))

  # Marginal likelihood of p1968
  k_like <- grid %>%
    group_by(k) %>%
    summarize(marg_like = sum(norm_posterior))

  return(list(grid = grid, r_like = r_like,k_like = k_like))
}
