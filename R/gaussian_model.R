gaussian_model <- function(dimension){
  
  # construct standard Gaussian prior distribution
  prior_mean <- rep(0, dimension)
  prior_cov <- diag(dimension)
  logprior <- function(x){
    return(mvnpdf(x, prior_mean, prior_cov))  
  }
  gradlogprior <- function(x){
    return(prior_mean - x)
  }
  sampleprior <- function(nparticles){
    return(mvnrnd(nparticles, prior_mean, prior_cov))
  }
  prior <- list(logdensity = logprior, gradlogdensity = gradlogprior, rinit = sampleprior)
  
  # construct Gaussian log-likelihood
  like_obs <- rep(14.25, dimension)
  like_rho <- 0.5
  like_cov <- like_rho * matrix(1, dimension, dimension) + diag(1 - like_rho, dimension)
  like_rooti <- t(solve(chol(like_cov)))
  like_rootisum <- sum(log(diag(like_rooti)))
  loglikelihood <- function(x){
    return(mvnpdf_chol(x, like_obs, like_rooti, like_rootisum))
  }
  like_precision <- solve(like_cov)
  gradloglikelihood <- function(x){
    return((like_obs-x) %*% like_precision)
  }
  likelihood <- list(logdensity = loglikelihood, gradlogdensity = gradloglikelihood)
  
  gaussian_model <- list(prior = prior, likelihood = likelihood)
  
  return(gaussian_model)
}
