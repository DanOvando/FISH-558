pella.tom.gridmodel <- function(g,grid,z = 2.39,dat){

  r <- grid$r[g]

  k <- grid$k[g]

  tau <- grid$tau[g]

  dat$nhat[1] <- k

  dat$sigma <- sqrt(dat$u^2 + tau^2)

  for (t in 2:dim(dat)[2])
  {
    last.n <- dat$nhat[t-1]

    dat$nhat[t] <- pmax(last.n + (last.n * r) * (1- (last.n/k)^z) - dat$catch[t-1],0.0001)
  }

  dat$likelihood <- 1/dat$sigma * exp( - (log(dat$nhat) - log (dat$n))^2/(2*dat$sigma^2))

  sum(-log(1/dat$sigma) +  (log(dat$nhat) - log (dat$n))^2/(2*dat$sigma^2),na.rm = T)

  output <- data.frame(grid[g,],likelihood = prod(dat$likelihood, na.rm = T),nll = -sum(log(dat$likelihood),na.rm = T))

  return(output)
}
