% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_gibbsflow_ais.R
\name{run_gibbsflow_ais}
\alias{run_gibbsflow_ais}
\title{Run Gibbs flow annealed importance sampler}
\usage{
run_gibbsflow_ais(prior, likelihood, nparticles, lambda, timegrid, mcmc)
}
\arguments{
\item{prior}{list with keys: 
\code{logdensity} evaluates log prior density, 
\code{gradlogdensity} returns its gradient, 
\code{rinit} samples from the prior distribution}

\item{likelihood}{list with keys: 
\code{logdensity} samples from proposal, 
\code{gradlogdensity} returns its gradient}

\item{nparticles}{number of particles}

\item{lambda}{vector describing tempering schedule}

\item{timegrid}{vector describing numerical integration times}

\item{mcmc}{list with keys: 
\code{choice} specifies type of MCMC method, 
\code{parameters} specifies algorithmic tuning parameters,
\code{nmoves} specifies number of MCMC move per temperature}
}
\value{
list with keys: 
\code{xtrajectory} trajectories, 
\code{xparticles} particles at terminal time,
\code{ess} effective sample size, 
\code{log_normconst} log normalizing constant,
\code{acceptprob} MCMC acceptance probabilities
}
\seealso{
\code{\link{run_gibbsflow_smc}} if resampling is desired
}
