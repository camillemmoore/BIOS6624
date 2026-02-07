###########################################################
## Program Name: 3_FullModel_Multivariable_Worksheet4.R  ##
## Purpose: Multivariable model for PSA                  ##
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


####Run a linear model just to orient us using standard methods
fullmodel <-lm(lpsa ~ cavol + wt + age + bph + cappen + svi + grade7 + grade8 , data=psaclean2)

sink("Output/WS4_FullModelOutput.txt")
print(fullmodel)
summary(fullmodel)
confint(fullmodel)
sink()



####Bayesian Analysis
# Compile the Stan program for linear regression
mod <- cmdstan_model('Code/STAN/linear_regression_half_normal.stan')

# Outcome data
y <- psaclean2$lpsa

# First fit full model
  mod_name <- 'full'
  X <- model.matrix(~ cavol + wt + age + bph + svi + cappen + grade7 + grade8, data = psaclean2) 
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
  res_full <- summary_table[!(summary_table$Parameter %in% c('beta[1]', 'sigma')),]
  res_full$param_name <- colnames(X)[-1]
  
# Full model diagnostic plots
  draws <- as_draws_array(fit$draws())  # dimensions: iterations x chains x parameters
  draws_df <- as_draws_df(fit$draws())  # tidy format for bayesplot / ggplot
  params <- c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "beta[6]", "beta[7]", "beta[8]", "beta[9]","sigma")
  
  pdf('Output/WS4_DiagnosticPlots.pdf')
  # Trace plots
  # Trace plots show whether chains are mixing well and exploring the parameter space.
  # Chains should mix well (lines overlap)
  # No long trends (no “stickiness”)
  # Ideally all chains overlap around the same mean
  
  mcmc_trace(draws, pars = params)
  
  # Posterior Density Plots
  # Overlaid densities from multiple chains should match closely
  # If densities differ substantially between chains → poor mixing
  mcmc_dens_overlay(draws, pars = params)
  
  # Auto-correlation 
  # Autocorrelation should decay quickly
  # High autocorrelation → may need more iterations or thinning
  mcmc_acf(draws, pars = params)
  
  
  dev.off()
  
  # Diagnostics from cmdstan
  fit$cmdstan_diagnose()
  
  # R-hat and ESS diagnostics
  # R-hat ≈ 1.00 → converged
  # ESS should be reasonably large (e.g., >100–200 per parameter)
  summary_stats <- summarise_draws(draws)
  
  summary_stats %>%
    filter(variable %in% params) %>%
    dplyr::select(variable, rhat, ess_bulk, ess_tail)
  
  
  
##### function to generate info for WS4
##### save posterior mean, 95%HPDI, DIC, WAIC, LOO-CV for when each var is removed
exclude_fits <- function(y, X, var_name){
  
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
  
  return(data.frame(
    DIC = DIC,
    WAIC=waic_res$estimates['waic', 'Estimate'],
    LOOCV=loo_res$estimates['looic', 'Estimate'],
    variable=var_name
    ))
  
}

res_exclude <- NULL
# cavol 
X <- model.matrix(~ wt + age + bph + svi + cappen + grade7 + grade8, data = psaclean2) 

res_exclude <- rbind(res_exclude, 
                   exclude_fits(y, X, 'cavol removed'))

# wt 
X <- model.matrix(~ cavol + age + bph + svi + cappen + grade7 + grade8, data = psaclean2) 

res_exclude <- rbind(res_exclude, 
                   exclude_fits(y, X, 'wt removed'))

# age 
X <- model.matrix(~ cavol + wt + bph + svi + cappen + grade7 + grade8, data = psaclean2) 

res_exclude <- rbind(res_exclude, 
                   exclude_fits(y, X, 'age removed'))


# bph
X <- model.matrix(~ cavol + wt + age + svi + cappen + grade7 + grade8, data = psaclean2) 

res_exclude <- rbind(res_exclude, 
                   exclude_fits(y, X, 'bph removed'))


# svi 
X <- model.matrix(~ cavol + wt + age + bph + cappen + grade7 + grade8, data = psaclean2) 

res_exclude <- rbind(res_exclude, 
                   exclude_fits(y, X, 'svi removed'))


# cappen 
X <- model.matrix(~ cavol + wt + age + bph + svi + grade7 + grade8, data = psaclean2) 

res_exclude <- rbind(res_exclude, 
                   exclude_fits(y, X, 'cappen removed'))

# gleason
X <- model.matrix(~ cavol + wt + age + bph + svi + cappen, data = psaclean2) 

res_exclude <- rbind(res_exclude, 
                   exclude_fits(y, X, 'gleason removed'))

# Save Results
write.xlsx(list(res_full, res_exclude), './Output/WS4_multivariable_fit_summary.xlsx')


