#' Process whale SIR
#'
#' \code{process.whale.SIR} takes the output of the
#' whale SIR and makes plots and diagnostics
#' @param fitted.whales dataframe of SIR results
#' @param dat the raw data
#' @param catch.dat the catch data
#' @param rungroup name of run to append to figures
#' @param runfolder the folder to store results and figures
#' @export

process.whale.SIR <- function(fitted.whales, dat, catch.dat, rungroup, runfolder, use.catch = T)
{

  posterior.distributions <- fitted.whales %>%
                                rename(s.a = s.rest) %>%
                                gather('variable','value',s.0:f.max) %>%
    group_by(variable) %>%
    mutate(rng = diff(range(value)), bnwdth = rng/20) %>%
    ungroup()

  posterior.distributions.plot  <- posterior.distributions %>%
    ggplot(aes(value, fill = variable)) +
    scale_fill_brewer(guide = F, palette = 'Spectral') +
    geom_histogram(aes(y = ..density..) ,alpha = 0.6) +
    geom_density(alpha = 0.6) +
    facet_wrap(~variable, scales = 'free') +
    theme_bw()

  ggsave(paste(runfolder,'/',rungroup,'-posterior distributions.pdf', sep = ''), plot = posterior.distributions.plot, height = 6, width = 6)

  check.SIR <- fitted.whales %>%
    ungroup() %>%
    mutate(param.space = paste(s.0,s.rest,K,f.max,sep = '-')) %>%
    group_by(param.space) %>%
    summarize(num.samples = length(s.0)) %>%
    ungroup() %>%
    mutate(perc.samples = num.samples / sum(num.samples))

  check.SIR.plot <- check.SIR %>%
    ggplot(aes(perc.samples)) +
    geom_histogram() +
    geom_vline(aes(xintercept = 0.05)) +
    scale_x_continuous(labels = percent) +
    ylab('# of Unique Parameter Vectors')

  ggsave(paste(runfolder,'/',rungroup,'-SIR Diagnostic.pdf', sep = ''), plot = check.SIR.plot, height = 5, width = 5)

  #   whale.fit <- make.whales(dat = dat, catch.dat = catch.dat, s.0 = mean(fitted.whales$s.0),
  #                            s.rest = mean(fitted.whales$s.rest), f.max = mean(fitted.whales$f.max),
  #                            K = mean(fitted.whales$K), use.catch = use.catch)

  simmed.whales <- sim.whales(dat = dat, catch.dat = catch.dat,
                              possible.whales = fitted.whales, use.catch = use.catch)

  labels <- seq(min(dat$year), max(dat$year),  length.out = 10)

  years <- unique(dat$year)

  breaks <- seq(min(years),max(years), length.out = 10)

  test <- as.factor(simmed.whales$year)

  posterior.predicted.whales <- simmed.whales %>%
    ggplot(aes(factor(year),predicted.whales, fill = 'Predicted Whales')) +
    geom_violin(alpha = 0.5) +
    geom_point(aes(factor(year),(abundance), fill = 'Observed Whales'), shape=  21,size = 3) +
    scale_fill_discrete(name = element_blank()) +
    #   scale_x_date(labels = date_format("%Y")) +
    scale_x_discrete(labels = labels, breaks = breaks) +
    theme_bw() +
    theme(legend.position = 'top') +
    xlab('Year') +
    ylab('Whales')

  ggsave(paste(runfolder,'/',rungroup,'- Distribution of predicted whales.pdf', sep = ''), plot = posterior.predicted.whales, height = 5, width = 8)

  summarized.whales <- simmed.whales %>%
    group_by(year) %>%
    summarize(lower.5 = quantile(predicted.whales,probs = 0.05, na.rm = T),
              upper.95 =quantile(predicted.whales,probs = 0.95, na.rm = T),
              median.whales = median(predicted.whales)) %>%
    left_join(dat, by = 'year')

  whale.summary.plot <- summarized.whales %>%
    ggplot(aes(x = year)) +
    geom_ribbon(aes(ymin = lower.5, ymax = upper.95, fill = '90% Range'), alpha = 0.25) +
    geom_line(aes(year, median.whales, color = 'Median'), linetype = 'longdash') +
    geom_point(aes(year, abundance, color = 'Observed'), shape = 19) +
    scale_color_manual(name = element_blank(), values = c('lightseagreen','tomato') ) +
    scale_fill_manual(name = element_blank(), values = 'grey50') +
    xlab('Year') +
    ylab('Numbers of Whales')+
    theme_bw()

  ggsave(paste(runfolder,'/',rungroup,'- Median distribution of predicted whales.pdf', sep = ''), plot = whale.summary.plot, height = 5, width = 8)


  return(list(posterior.distributions = posterior.distributions,posterior.distributions.plot = posterior.distributions.plot, simmed.whales = simmed.whales,
              summarized.whales = summarized.whales,check.SIR = check.SIR, whale.summary.plot = whale.summary.plot,
              posterior.predicted.whales.plot = posterior.predicted.whales, check.SIR.plot = check.SIR.plot))
}