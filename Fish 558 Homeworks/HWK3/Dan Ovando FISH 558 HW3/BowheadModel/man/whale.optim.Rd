% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/whale_optim.R
\name{whale.optim}
\alias{whale.optim}
\title{Do a simple SIR to test validity of parameters}
\usage{
whale.optim(pars, dat, catch.dat)
}
\arguments{
\item{dat}{the raw whale data}

\item{s.0.}{}

\item{s.rest}{adult survival}

\item{K}{carrying capacity}

\item{f.max}{max fecundity}
}
\description{
\code{whale.optim} returns a posterior of calf and adult
survival, max fecundity, and K, such that the population
does not crash given an observed catch history
}

