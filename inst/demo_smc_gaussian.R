rm(list = ls())
library(tictoc)
library(ggplot2)

# model
dimension <- 10
model <- gaussian_model(dimension)

# SMC settings
nparticles <- 2^10
nsteps <- 20
lambda <- seq(0, 1, length.out = nsteps)
mcmc <- list()
mcmc$choice <- "rwmh"
mcmc$parameters$stepsize <- 0.1
mcmc$parameters$nsteps <- 10
mcmc$nmoves <- 5

# run SMC
tic("SMC runtime")
  smc <- run_smc(model$prior, model$likelihood, nparticles, lambda, mcmc)
toc()

# ess plot
ess.df <- data.frame(time = 1:nsteps, ess = smc$ess * 100 / nparticles)
ggplot(ess.df, aes(x = time, y = ess)) + geom_point()

# normalizing constant plot
normconst.df <- data.frame(time = 1:nsteps, normconst = smc$log_normconst)
ggplot(normconst.df, aes(x = time, y = normconst)) + geom_point()
smc$log_normconst[nsteps]

