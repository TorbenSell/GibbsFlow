rm(list = ls())
library(GibbsFlow)
library(tictoc)
library(ggplot2)

# prior
prior <- list()
prior$logdensity <- function(x) as.numeric(gaussian_logprior(x))
prior$gradlogdensity <- function(x) gaussian_gradlogprior(x)
prior$rinit <- function(n) gaussian_sampleprior(n) 

# likelihood
likelihood <- list()
likelihood$logdensity <- function(x) as.numeric(gaussian_loglikelihood(x))
likelihood$gradlogdensity <- function(x) gaussian_gradloglikelihood(x)
compute_gibbsflow <- function(stepsize, lambda, derivative_lambda, x, logdensity) gaussian_gibbsflow(stepsize, lambda, derivative_lambda, x, logdensity)

# posterior
posterior <- list()
posterior$mean <- function(lambda) gaussian_posterior_mean(lambda)
posterior$cov <- function(lambda) gaussian_posterior_cov(lambda)
posterior$log_normconst <- function(lambda) gaussian_log_normconst(lambda)
posterior$logdensity <- function(x, lambda) gaussian_logposterior(x, lambda)

# smc settings
nparticles <- 2^10
nsteps <- 20
timegrid <- seq(0, 1, length.out = nsteps)
exponent <- 1
lambda <- timegrid^exponent
derivative_lambda <- exponent * timegrid^(exponent - 1)

# run sampler
set.seed(7)
tic()
  smc <- run_gibbsflow_sis(prior, likelihood, nparticles, timegrid, lambda, derivative_lambda, compute_gibbsflow)
toc()

# ess plot
ess.df <- data.frame(time = 1:nsteps, ess = smc$ess * 100 / nparticles)
ggplot(ess.df, aes(x = time, y = ess)) + geom_line() + 
  labs(x = "time", y = "ESS%") + ylim(c(0, 100))

# normalizing constant plot
normconst.df <- data.frame(time = 1:nsteps, normconst = smc$log_normconst)
true_normconst.df <- data.frame(time = 1:nsteps, normconst = sapply(lambda, posterior$log_normconst))
ggplot() + geom_line(data = normconst.df, aes(x = time, y = normconst), colour = "blue") + 
  geom_line(data = true_normconst.df, aes(x = time, y = normconst), colour = "red") +
  labs(x = "time", y = "log normalizing constant")

# plot density function and trajectory
xgrid <- seq(-3, 7, length.out = 200)
ygrid <- seq(-3, 7, length.out = 200)
grid <- expand.grid(xgrid = xgrid, ygrid = ygrid)
prior_values <- rep(0, nrow(grid))
posterior_values <- rep(0, nrow(grid))
for (igrid in 1:nrow(grid)){
  prior_values[igrid] <- exp(prior$logdensity(as.matrix(grid[igrid, ], nrow = 1)))
  posterior_values[igrid] <- exp(posterior$logdensity(as.matrix(grid[igrid, ], nrow = 1), 1))
}
plot_prior <- cbind(grid, prior_values)
plot_posterior <- cbind(grid, posterior_values)
plot_particles <- data.frame(x = smc$xparticles[, 1], y = smc$xparticles[, 2])
ggplot() + stat_contour(data = plot_prior, aes(x = xgrid, y = ygrid, z = prior_values), colour = "red") +
  stat_contour(data = plot_posterior, aes(x = xgrid, y = ygrid, z = posterior_values)) +
  geom_point(data = plot_particles, aes(x = x, y = y), colour = "black")


