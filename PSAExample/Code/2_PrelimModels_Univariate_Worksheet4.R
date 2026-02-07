###########################################################
## Program Name: 2_PrelimModels_Univariate_Worksheet3.R  ##
## Purpose: Investigate independent effects of PSA for   ##
##           each variable separately.                   ## 
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


####Bayesian Analysis
# Compile the Stan program for linear regression
mod <- cmdstan_model('Code/STAN/linear_regression_half_normal.stan')

# Outcome data
y <- psaclean2$lpsa

# First intercept only
mod_name <- 'Intercept only'
X <- model.matrix(~ 1, data = psaclean2) 
#These are variables needed in Stan
#N is the number of observations in your dataset and P is the number of columns in the design matrix
N <- nrow(X)
P <- ncol(X)

## Hyperparameters for prior distributions. These change depending on your assumptions
## For this problem the model error SD (sigma) is a half normal and the coefficients
## are normal distribution with a mean and SD, which can be written flexibly 
## with a mean vector (m below) and a SD vector (s below).  

m <- c(1, rep(0, (P-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
s <- rep(100,P) # SD in the prior on regression coefficients --> variance 100^2
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
summary_table$variable <- mod_name
res_final <- summary_table[(summary_table$Parameter %in% c('beta[1]')),]


##### function to generate info for WS4
##### save posterior mean, 95%HPDI, DIC, WAIC, LOO-CV for each variable
univariate_fits <- function(y, X, var_name){
  
  #These are variables needed in Stan
  #N is the number of observations in your dataset and P is the number of columns in the design matrix
  N <- nrow(X)
  P <- ncol(X)
  
  ## Hyperparameters for prior distributions. These change depending on your assumptions
  ## For this problem the model error SD (sigma) is a half normal and the coefficients
  ## are normal distribution with a mean and SD, which can be written flexibly 
  ## with a mean vector (m below) and a SD vector (s below).  
  
  m <- c(1, rep(0, (P-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
  s <- rep(100,P) # SD in the prior on regression coefficients --> variance 100^2
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
  summary_table$variable <- var_name
  
  return(summary_table[!(summary_table$Parameter %in% c('beta[1]', 'sigma')),])
  
}

# cavol 
X <- model.matrix(~ cavol, data = psaclean2) 

res_final <- rbind(res_final, 
                     univariate_fits(y, X, 'cavol'))

# wt 
X <- model.matrix(~ wt, data = psaclean2) 

res_final <- rbind(res_final, 
                   univariate_fits(y, X, 'wt'))

# age 
X <- model.matrix(~ age, data = psaclean2) 

res_final <- rbind(res_final, 
                   univariate_fits(y, X, 'age'))


# bph
X <- model.matrix(~ bph, data = psaclean2) 

res_final <- rbind(res_final, 
                   univariate_fits(y, X, 'bph'))


# svi 
X <- model.matrix(~ svi, data = psaclean2) 

res_final <- rbind(res_final, 
                   univariate_fits(y, X, 'svi'))


# cappen 
X <- model.matrix(~ cappen, data = psaclean2) 

res_final <- rbind(res_final, 
                   univariate_fits(y, X, 'cappen'))

# gleason
X <- model.matrix(~ grade7+grade8, data = psaclean2) 

res_final <- rbind(res_final, 
                   univariate_fits(y, X, 'gleason'))

# Save Results
write.xlsx(res_final, './Output/WS4_univariate_fit_summary.xlsx')




