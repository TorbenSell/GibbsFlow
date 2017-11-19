#' @rdname run_gibbsflow_sis 
#' @title Run Gibbs flow sequential importance sampler
#' @param prior list with keys: 
#' \code{logdensity} evaluates log prior density, 
#' \code{gradlogdensity} returns its gradient, 
#' \code{rinit} samples from the prior distribution
#' @param likelihood list with keys: 
#' \code{logdensity} samples from proposal, 
#' \code{gradlogdensity} returns its gradient
#' @param nparticles number of particles
#' @param lambda vector describing tempering schedule
#' @param timegrid vector describing numerical integration times 
#' @return list with keys: 
#' \code{xtrajectory} trajectories, 
#' \code{xparticles} particles at terminal time,
#' \code{ess} effective sample size, 
#' \code{log_normconst} log normalizing constant, 
#' @seealso \code{\link{run_gibbsflow_smc}} if resampling is desired
#' @export
run_gibbsflow_sis <- function(prior, likelihood, nparticles, lambda, timegrid, compute_gibbsflow){
  # initialization
  xparticles <- prior$rinit(nparticles)
  previous_logdensity <- prior$logdensity(xparticles)
  
  # pre-allocate
  dimension <- ncol(xparticles)
  nsteps <- length(lambda) # same length as timegrid
  stepsize <- diff(timegrid) 
  xtrajectory <- array(dim = c(nparticles, dimension, nsteps))
  xtrajectory[ , , 1] <- xparticles
  logweights <- rep(0, nparticles)
  ess <- rep(0, nsteps)
  ess[1] <- nparticles
  log_normconst <- rep(0, nsteps)
  
  for (istep in 2:nsteps){
    # gibbs flow move
    output_flow <- compute_gibbsflow(stepsize[istep-1], timegrid[istep-1], xparticles, previous_logdensity)
    xparticles <- output_flow$xparticles 
    jacobian_det <- output_flow$jacobian_det 

    # weight 
    current_logdensity <- prior$logdensity(xparticles) + lambda[istep] * likelihood$logdensity(xparticles)
    logweights <- logweights + current_logdensity - previous_logdensity + log(abs(jacobian_det))
    maxlogweights <- max(logweights)
    weights <- exp(logweights - maxlogweights)
    normweights <- weights / sum(weights)  
    
    # compute effective sample size
    ess[istep] <- 1 / sum(normweights^2)
    
    # compute normalizing constant
    log_normconst[istep] <- log(mean(weights)) + maxlogweights
    
    # store trajectory
    xtrajectory[ , , istep] <- xparticles
    previous_logdensity <- current_logdensity
    
  }
  
  return(list(xtrajectory = xtrajectory, xparticles = xparticles, ess = ess, 
              log_normconst = log_normconst))
}

