###########################################################
## Program Name: 5_SVI_Interaction_Models_Worksheet5.R   ##
## Purpose: Investigate SVI Interactions with Wt and Cavol##
## Created by: Camille Moore                             ##
###########################################################

# Dependencies
library(cmdstanr)
library(bayesplot)
library(posterior)
library(bayestestR)
library(mcmcse)
library(loo)
library(dplyr)
library(tibble)

setwd("/Users/mooreca/Documents/BIOS6624/PSAExample/")
psaclean2 <- read.csv("DataProcessed/psaclean2.csv")


####Run linear models just to orient us using standard methods
cavol_model <-lm(lpsa ~ cavol + svi + cavolsviint, data=psaclean2)
weight_model <-lm(lpsa ~ wt + svi + wtsviint, data=psaclean2)

sink("Output/WS5_SVI_Interaction_Output.txt")
print(cavol_model)
summary.lm(cavol_model)
confint(cavol_model)

print(weight_model)
summary.lm(weight_model)
confint(weight_model)
sink()

####Bayesian Analysis
# Compile the Stan program for linear regression
mod <- cmdstan_model('Code/STAN/linear_regression_half_normal.stan')

# Outcome data
y <- psaclean2$lpsa


# First fit full model for Cavol
mod_name <- 'cavol_svi_interaction'
X <- model.matrix(~ cavol + svi + cavolsviint, data = psaclean2) 
#These are variables needed in Stan
#N is the number of observations in your dataset and P is the number of columns in the design matrix
N <- nrow(X)
P <- ncol(X)

## Hyperparameters for prior distributions. These change depending on your assumptions
## For this problem the model error SD (sigma) is a half normal and the coefficients
## are normal distribution with a mean and SD, which can be written flexibly 
## with a mean vector (m below) and a SD vector (s below).  

m <- c(1, rep(0, (P-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
s <- rep(sqrt(1000),P) # SD in the prior on regression coefficients --> variance 100^2
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
  iter_warmup = 2000,
  iter_sampling = 10000,
  seed = 123
)

# Table of Results of Interest from the Posterior
# Extract posterior draws from cmdstanr fit
draws <- fit$draws()  # cmdstanr fit

draws_mat <- as_draws_matrix(fit$draws())
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

####Model Fit Stats
# Get LOO-CV and WAIC
loglik_mat <- as_draws_matrix(fit$draws("log_lik"))  # iterations x N

loo_res <- loo(loglik_mat)
waic_res <- waic(loglik_mat)

# Compute DIC by hand
# Deviance per iteration: -2 * log-likelihood summed over observations
D_theta <- -2 * rowSums(loglik_mat)
mean_D <- mean(D_theta)

# Deviance at posterior mean
beta_mean <- colMeans(draws_mat[, grep("^beta\\[", colnames(draws_mat))])
sigma_mean <- mean(draws_mat[, "sigma"])
D_bar_theta <- -2 * sum(dnorm(y, mean = X %*% beta_mean, sd = sigma_mean, log = TRUE))

DIC <- 2 * mean_D - D_bar_theta

summary_table$DIC <- DIC
summary_table$WAIC <- waic_res$estimates['waic', 'Estimate']
summary_table$LOOCV <- loo_res$estimates['looic', 'Estimate']
summary_table$model <- mod_name
summary_table$variable <- c(colnames(X), 'sigma')
cavol_int_result <- summary_table


### cavol svi model without interaction
mod_name <- 'cavol_svi_noint'
X <- model.matrix(~ cavol + svi, data = psaclean2) 
#These are variables needed in Stan
#N is the number of observations in your dataset and P is the number of columns in the design matrix
N <- nrow(X)
P <- ncol(X)

## Hyperparameters for prior distributions. These change depending on your assumptions
## For this problem the model error SD (sigma) is a half normal and the coefficients
## are normal distribution with a mean and SD, which can be written flexibly 
## with a mean vector (m below) and a SD vector (s below).  

m <- c(1, rep(0, (P-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
s <- rep(sqrt(1000),P) # SD in the prior on regression coefficients --> variance 100^2
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
  iter_warmup = 2000,
  iter_sampling = 10000,
  seed = 123
)

# Table of Results of Interest from the Posterior
# Extract posterior draws from cmdstanr fit
draws <- fit$draws()  # cmdstanr fit

draws_mat <- as_draws_matrix(fit$draws())
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

####Model Fit Stats
# Get LOO-CV and WAIC
loglik_mat <- as_draws_matrix(fit$draws("log_lik"))  # iterations x N

loo_res <- loo(loglik_mat)
waic_res <- waic(loglik_mat)

# Compute DIC by hand
# Deviance per iteration: -2 * log-likelihood summed over observations
D_theta <- -2 * rowSums(loglik_mat)
mean_D <- mean(D_theta)

# Deviance at posterior mean
beta_mean <- colMeans(draws_mat[, grep("^beta\\[", colnames(draws_mat))])
sigma_mean <- mean(draws_mat[, "sigma"])
D_bar_theta <- -2 * sum(dnorm(y, mean = X %*% beta_mean, sd = sigma_mean, log = TRUE))

DIC <- 2 * mean_D - D_bar_theta

summary_table$DIC <- DIC
summary_table$WAIC <- waic_res$estimates['waic', 'Estimate']
summary_table$LOOCV <- loo_res$estimates['looic', 'Estimate']
summary_table$model <- mod_name
summary_table$variable <- c(colnames(X), 'sigma')
cavol_svi_reduced_result <- summary_table


### weight x svi analysis

# First fit full model for Cavol
mod_name <- 'weight_svi_interaction'
X <- model.matrix(~ wt + svi + wtsviint, data = psaclean2) 
#These are variables needed in Stan
#N is the number of observations in your dataset and P is the number of columns in the design matrix
N <- nrow(X)
P <- ncol(X)

## Hyperparameters for prior distributions. These change depending on your assumptions
## For this problem the model error SD (sigma) is a half normal and the coefficients
## are normal distribution with a mean and SD, which can be written flexibly 
## with a mean vector (m below) and a SD vector (s below).  

m <- c(1, rep(0, (P-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
s <- rep(sqrt(1000),P) # SD in the prior on regression coefficients --> variance 100^2
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
  iter_warmup = 2000,
  iter_sampling = 10000,
  seed = 123
)

# Table of Results of Interest from the Posterior
# Extract posterior draws from cmdstanr fit
draws <- fit$draws()  # cmdstanr fit

draws_mat <- as_draws_matrix(fit$draws())
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

####Model Fit Stats
# Get LOO-CV and WAIC
loglik_mat <- as_draws_matrix(fit$draws("log_lik"))  # iterations x N

loo_res <- loo(loglik_mat)
waic_res <- waic(loglik_mat)

# Compute DIC by hand
# Deviance per iteration: -2 * log-likelihood summed over observations
D_theta <- -2 * rowSums(loglik_mat)
mean_D <- mean(D_theta)

# Deviance at posterior mean
beta_mean <- colMeans(draws_mat[, grep("^beta\\[", colnames(draws_mat))])
sigma_mean <- mean(draws_mat[, "sigma"])
D_bar_theta <- -2 * sum(dnorm(y, mean = X %*% beta_mean, sd = sigma_mean, log = TRUE))

DIC <- 2 * mean_D - D_bar_theta

summary_table$DIC <- DIC
summary_table$WAIC <- waic_res$estimates['waic', 'Estimate']
summary_table$LOOCV <- loo_res$estimates['looic', 'Estimate']
summary_table$model <- mod_name
summary_table$variable <- c(colnames(X), 'sigma')
weight_int_result <- summary_table


### wt svi model without interaction
mod_name <- 'wt_svi_noint'
X <- model.matrix(~ wt + svi, data = psaclean2) 
#These are variables needed in Stan
#N is the number of observations in your dataset and P is the number of columns in the design matrix
N <- nrow(X)
P <- ncol(X)

## Hyperparameters for prior distributions. These change depending on your assumptions
## For this problem the model error SD (sigma) is a half normal and the coefficients
## are normal distribution with a mean and SD, which can be written flexibly 
## with a mean vector (m below) and a SD vector (s below).  

m <- c(1, rep(0, (P-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
s <- rep(sqrt(1000),P) # SD in the prior on regression coefficients --> variance 100^2
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
  iter_warmup = 2000,
  iter_sampling = 10000,
  seed = 123
)

# Table of Results of Interest from the Posterior
# Extract posterior draws from cmdstanr fit
draws <- fit$draws()  # cmdstanr fit

draws_mat <- as_draws_matrix(fit$draws())
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

####Model Fit Stats
# Get LOO-CV and WAIC
loglik_mat <- as_draws_matrix(fit$draws("log_lik"))  # iterations x N

loo_res <- loo(loglik_mat)
waic_res <- waic(loglik_mat)

# Compute DIC by hand
# Deviance per iteration: -2 * log-likelihood summed over observations
D_theta <- -2 * rowSums(loglik_mat)
mean_D <- mean(D_theta)

# Deviance at posterior mean
beta_mean <- colMeans(draws_mat[, grep("^beta\\[", colnames(draws_mat))])
sigma_mean <- mean(draws_mat[, "sigma"])
D_bar_theta <- -2 * sum(dnorm(y, mean = X %*% beta_mean, sd = sigma_mean, log = TRUE))

DIC <- 2 * mean_D - D_bar_theta

summary_table$DIC <- DIC
summary_table$WAIC <- waic_res$estimates['waic', 'Estimate']
summary_table$LOOCV <- loo_res$estimates['looic', 'Estimate']
summary_table$model <- mod_name
summary_table$variable <- c(colnames(X), 'sigma')
weight_svi_reduced_result <- summary_table


library(openxlsx)

# Save Results
write.xlsx(list(cavol_int_result, cavol_svi_reduced_result, weight_int_result, weight_svi_reduced_result
), './Output/WS5_svi_interaction_model_summary.xlsx')


