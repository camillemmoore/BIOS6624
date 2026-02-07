############################################################
## Program Name: 2_PrelimModels_Priors_and_MCMC_Lecture.R ##
## Purpose: Investigate the effect of age on PSA          ##
##           using different priors                       ## 
## Created by: Camille Moore                              ##
############################################################

library(cmdstanr)
library(posterior)
library(bayestestR)
library(mcmcse)
library(loo)
library(dplyr)
library(tibble)

setwd("/Users/mooreca/Documents/BIOS6624/PSAExample/")
psaclean2 <- read.csv("DataProcessed/psaclean2.csv")

###########################################################################
#################LOG PSA AND AGE ANALYSIS##################################
###########################################################################

####Run a linear model just to orient us using standard methods
agemodel <-lm(lpsa ~ age, data=psaclean2)

sink("Output/Lecture_AgeLinearModelOutput.txt")
print(agemodel)
summary(agemodel)
confint(agemodel)
sink()

####Bayesian Analysis

# Compile the Stan program
# Once the stan file is written it can be reused multiple times
mod <- cmdstan_model('Code/STAN/linear_regression_half_normal.stan')


# Prepare your data to pass to Stan
# Outcome data
y <- psaclean2$lpsa

# Now we set up the design matrix in our linear regression
X <- model.matrix(~ age, data = psaclean2) 

#These are variables needed in Stan
#N is the number of observations in your dataset and P is the number of columns in the design matrix
N <- nrow(X)
P <- ncol(X)

## Hyperparameters for prior distributions. These change depending on your assumptions
## For this problem the model error SD (sigma) is a half normal and the coefficients
## are normal distribution with a mean and SD, which can be written flexibly 
## with a mean vector (m below) and a SD vector (s below).  

################################################################################
# Non-Informative Prior: N(0, 500^2)
################################################################################
m <- c(1, rep(0, (P-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
s <- rep(500,P) # SD in the prior on regression coefficients --> variance 250000
sigma_sd <- 100

# create data list to pass to STAN
data_list <- list(
  N = N,
  P = P,
  X = X,
  y = y,
  prior_mean = m,
  prior_sd = s,
  sigma_prior_sd = sigma_sd
)

# Fit model
fit <- mod$sample(
  data = data_list,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 25000,
  seed = 123
)

# 1. Extract posterior draws from cmdstanr fit
draws <- fit$draws()  # cmdstanr fit

draws_mat <- as_draws_matrix(fit$draws())
draws_df <- as_draws_df(fit$draws())  # tidy format for bayesplot / ggplot
params <- colnames(draws_mat)
# Exclude non-parameter columns: lp__ and log_lik[n]
params <- params[!grepl("lp__|log_lik", params)]


# Build summary table
summary_table <- lapply(params, function(p) {
  vals <- as.numeric(draws_mat[, p])
  
  # Monte Carlo standard error
  mcse_val <- mcmcse::mcse(vals)$se
  
  # Effective sample size
  ess_val <- ess_bulk(vals)
  
  # 95% HPDI
  hpd <- hdi(vals, ci = 0.95)
  
  tibble(
    Parameter = p,
    Estimate  = mean(vals),
    MCSE      = mcse_val,
    Std_Dev   = sd(vals),
    HPDI_2.5  = hpd$CI_low,
    HPDI_97.5 = hpd$CI_high,
    ESS       = ess_val
  )
}) %>% bind_rows()

print(summary_table)

# Plot the Posterior and make diagnostic plots
p1 <- ggplot(draws_df, aes(x=`beta[1]`))+geom_density()+theme_classic()+
  xlab('Intercept')+ylab('Density')

p2 <- ggplot(draws_df, aes(x=`beta[2]`))+geom_density()+theme_classic()+
  xlab('Age')+ylab('Density')

p3 <- ggplot(draws_df, aes(x=sigma))+geom_density()+theme_classic()+
  xlab('Sigma')+ylab('Density')

ggarrange(p1, p2, p3, nrow=3, ncol=1)


# Trace plots
# Trace plots show whether chains are mixing well and exploring the parameter space.
# Chains should mix well (lines overlap)
# No long trends (no “stickiness”)
# Ideally all chains overlap around the same mean
mcmc_trace(draws, pars = params)

# R-hat and ESS diagnostics
# R-hat ≈ 1.00 → converged
# ESS should be reasonably large (e.g., >100–200 per parameter)
fit$summary(variables=c("beta[1]", "beta[2]", "sigma"))

# Posterior Density Plots
# Overlaid densities from multiple chains should match closely
# If densities differ substantially between chains → poor mixing
mcmc_dens_overlay(draws, pars = params)

# Auto-correlation 
# Autocorrelation should decay quickly
# High autocorrelation → may need more iterations or thinning
mcmc_acf(draws, pars = params)

# Diagnostics from cmdstan
fit$cmdstan_diagnose()


### Evaluating the effect of age
# Posterior P for a 10% change in PSA with 1 year change in age 
prop.table(table(draws_df$`beta[2]`>log(1.1)|draws_df$`beta[2]`<log(0.9)))

# Posterior P for a 5% change in PSA with 1 year change in age 
prop.table(table(draws_df$`beta[2]`>log(1.05)|draws_df$`beta[2]`<log(0.95)))

# Get model fit statistics
loglik_mat <- as_draws_matrix(fit$draws("log_lik"))  # iterations x N

# LOO-IC
loo_res <- loo(loglik_mat)
print(loo_res)

# WAIC
waic_res <- waic(loglik_mat)
print(waic_res)

# Compute DIC by hand
# Deviance per iteration: -2 * log-likelihood summed over observations
D_theta <- -2 * rowSums(loglik_mat)
mean_D <- mean(D_theta)

# Deviance at posterior mean
beta_mean <- colMeans(draws_mat[, grep("^beta\\[", colnames(draws_mat))])
sigma_mean <- mean(draws_mat[, "sigma"])
D_bar_theta <- -2 * sum(dnorm(y, mean = X %*% beta_mean, sd = sigma_mean, log = TRUE))

DIC <- 2 * mean_D - D_bar_theta

print(DIC)

################################################################################
# Vague Prior: N(0, 100^2)
################################################################################
m <- c(1, rep(0, (P-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
s <- rep(100,P) # SD in the prior on regression coefficients --> variance 10000
sigma_sd <- 100

# create data list to pass to STAN
data_list <- list(
  N = N,
  P = P,
  X = X,
  y = y,
  prior_mean = m,
  prior_sd = s,
  sigma_prior_sd = sigma_sd
)

# Fit model
fit <- mod$sample(
  data = data_list,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 20000,
  seed = 123
)

# 1. Extract posterior draws from cmdstanr fit
draws <- fit$draws()  # cmdstanr fit

draws_mat <- as_draws_matrix(fit$draws())
draws_df <- as_draws_df(fit$draws())  # tidy format for bayesplot / ggplot

params <- colnames(draws_mat)
# Exclude non-parameter columns: lp__ and log_lik[n]
params <- params[!grepl("lp__|log_lik", params)]


# Build summary table
summary_table <- lapply(params, function(p) {
  vals <- as.numeric(draws_mat[, p])
  
  # Monte Carlo standard error
  mcse_val <- mcmcse::mcse(vals)$se
  
  # Effective sample size
  ess_val <- ess_bulk(vals)
  
  # 95% HPDI
  hpd <- hdi(vals, ci = 0.95)
  
  tibble(
    Parameter = p,
    Estimate  = mean(vals),
    MCSE      = mcse_val,
    Std_Dev   = sd(vals),
    HPDI_2.5  = hpd$CI_low,
    HPDI_97.5 = hpd$CI_high,
    ESS       = ess_val
  )
}) %>% bind_rows()

print(summary_table)

# Plot the Posterior
p1 <- ggplot(draws_df, aes(x=`beta[1]`))+geom_density()+theme_classic()+
  xlab('Intercept')+ylab('Density')

p2 <- ggplot(draws_df, aes(x=`beta[2]`))+geom_density()+theme_classic()+
  xlab('Age')+ylab('Density')

p3 <- ggplot(draws_df, aes(x=sigma))+geom_density()+theme_classic()+
  xlab('Sigma')+ylab('Density')

ggarrange(p1, p2, p3, nrow=3, ncol=1)

################################################################################
# Skeptical Prior: N(0, 10^2)
################################################################################
m <- c(1, rep(0, (P-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
s <- rep(10,P) # SD in the prior on regression coefficients --> variance 100
sigma_sd <- 100

# create data list to pass to STAN
data_list <- list(
  N = N,
  P = P,
  X = X,
  y = y,
  prior_mean = m,
  prior_sd = s,
  sigma_prior_sd = sigma_sd
)

# Fit model
fit <- mod$sample(
  data = data_list,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 25000,
  seed = 123
)

# 1. Extract posterior draws from cmdstanr fit
draws <- fit$draws()  # cmdstanr fit

draws_mat <- as_draws_matrix(fit$draws())
draws_df <- as_draws_df(fit$draws())  # tidy format for bayesplot / ggplot

params <- colnames(draws_mat)
# Exclude non-parameter columns: lp__ and log_lik[n]
params <- params[!grepl("lp__|log_lik", params)]


# Build summary table
summary_table <- lapply(params, function(p) {
  vals <- as.numeric(draws_mat[, p])
  
  # Monte Carlo standard error
  mcse_val <- mcmcse::mcse(vals)$se
  
  # Effective sample size
  ess_val <- ess_bulk(vals)
  
  # 95% HPDI
  hpd <- hdi(vals, ci = 0.95)
  
  tibble(
    Parameter = p,
    Estimate  = mean(vals),
    MCSE      = mcse_val,
    Std_Dev   = sd(vals),
    HPDI_2.5  = hpd$CI_low,
    HPDI_97.5 = hpd$CI_high,
    ESS       = ess_val
  )
}) %>% bind_rows()

print(summary_table)

# Plot the Posterior
p1 <- ggplot(draws_df, aes(x=`beta[1]`))+geom_density()+theme_classic()+
  xlab('Intercept')+ylab('Density')

p2 <- ggplot(draws_df, aes(x=`beta[2]`))+geom_density()+theme_classic()+
  xlab('Age')+ylab('Density')

p3 <- ggplot(draws_df, aes(x=sigma))+geom_density()+theme_classic()+
  xlab('Sigma')+ylab('Density')

ggarrange(p1, p2, p3, nrow=3, ncol=1)

################################################################################
# Optimistic Prior: N(10, 10^2)
################################################################################
m <- c(1, rep(10, (P-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
s <- rep(10,P) # SD in the prior on regression coefficients --> variance 100
sigma_sd <- 100

# create data list to pass to STAN
data_list <- list(
  N = N,
  P = P,
  X = X,
  y = y,
  prior_mean = m,
  prior_sd = s,
  sigma_prior_sd = sigma_sd
)

# Fit model
fit <- mod$sample(
  data = data_list,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 25000,
  seed = 123
)

# 1. Extract posterior draws from cmdstanr fit
draws <- fit$draws()  # cmdstanr fit

draws_mat <- as_draws_matrix(fit$draws())
draws_df <- as_draws_df(fit$draws())  # tidy format for bayesplot / ggplot

params <- colnames(draws_mat)
# Exclude non-parameter columns: lp__ and log_lik[n]
params <- params[!grepl("lp__|log_lik", params)]


# Build summary table
summary_table <- lapply(params, function(p) {
  vals <- as.numeric(draws_mat[, p])
  
  # Monte Carlo standard error
  mcse_val <- mcmcse::mcse(vals)$se
  
  # Effective sample size
  ess_val <- ess_bulk(vals)
  
  # 95% HPDI
  hpd <- hdi(vals, ci = 0.95)
  
  tibble(
    Parameter = p,
    Estimate  = mean(vals),
    MCSE      = mcse_val,
    Std_Dev   = sd(vals),
    HPDI_2.5  = hpd$CI_low,
    HPDI_97.5 = hpd$CI_high,
    ESS       = ess_val
  )
}) %>% bind_rows()

print(summary_table)

# Plot the Posterior
p1 <- ggplot(draws_df, aes(x=`beta[1]`))+geom_density()+theme_classic()+
  xlab('Intercept')+ylab('Density')

p2 <- ggplot(draws_df, aes(x=`beta[2]`))+geom_density()+theme_classic()+
  xlab('Age')+ylab('Density')

p3 <- ggplot(draws_df, aes(x=sigma))+geom_density()+theme_classic()+
  xlab('Sigma')+ylab('Density')

ggarrange(p1, p2, p3, nrow=3, ncol=1)

################################################################################
# Extreme Prior: N(10, 1^2)
################################################################################
m <- c(10, rep(10, (P-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
s <- rep(1,P) # SD in the prior on regression coefficients --> variance 1
sigma_sd <- 100

# create data list to pass to STAN
data_list <- list(
  N = N,
  P = P,
  X = X,
  y = y,
  prior_mean = m,
  prior_sd = s,
  sigma_prior_sd = sigma_sd
)

# Fit model
fit <- mod$sample(
  data = data_list,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 25000,
  seed = 123
)

# 1. Extract posterior draws from cmdstanr fit
draws <- fit$draws()  # cmdstanr fit

draws_mat <- as_draws_matrix(fit$draws())
draws_df <- as_draws_df(fit$draws())  # tidy format for bayesplot / ggplot

params <- colnames(draws_mat)
# Exclude non-parameter columns: lp__ and log_lik[n]
params <- params[!grepl("lp__|log_lik", params)]


# Build summary table
summary_table <- lapply(params, function(p) {
  vals <- as.numeric(draws_mat[, p])
  
  # Monte Carlo standard error
  mcse_val <- mcmcse::mcse(vals)$se
  
  # Effective sample size
  ess_val <- ess_bulk(vals)
  
  # 95% HPDI
  hpd <- hdi(vals, ci = 0.95)
  
  tibble(
    Parameter = p,
    Estimate  = mean(vals),
    MCSE      = mcse_val,
    Std_Dev   = sd(vals),
    HPDI_2.5  = hpd$CI_low,
    HPDI_97.5 = hpd$CI_high,
    ESS       = ess_val
  )
}) %>% bind_rows()

print(summary_table)

# Plot the Posterior
p1 <- ggplot(draws_df, aes(x=`beta[1]`))+geom_density()+theme_classic()+
  xlab('Intercept')+ylab('Density')

p2 <- ggplot(draws_df, aes(x=`beta[2]`))+geom_density()+theme_classic()+
  xlab('Age')+ylab('Density')

p3 <- ggplot(draws_df, aes(x=sigma))+geom_density()+theme_classic()+
  xlab('Sigma')+ylab('Density')

ggarrange(p1, p2, p3, nrow=3, ncol=1)
