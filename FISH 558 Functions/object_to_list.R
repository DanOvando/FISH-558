objects.to.list <- function(omit)
{

  files <- ls()

  files <- files[!(files %in% omit)]

  dat <- list()

  for (f in 1:length(files))
  {
    eval(parse(text = paste('dat$',files[f],'=',files[f], sep = '')))
  }

  return(dat)
}