############################################################################
# Multi-state model with capture-recapture analysis
# with Jolly-Seber model framework
# Code adapted from Kery and Schaub, BPA Book, Chapter 10
# Part 10.3.2

# Checking model diagnostics
# Same approach applies to all three models (additive, interaction 1, interaction 2)
# Example shown for the additive model

# Date: June 2026
# Authors: Dilsad Dagtekin and Dominik Behr
###########################################################################

## House keeping ----
##########################

rm(list = ls())
setwd("/path_to_files/")

## Load libraries ----
##########################

load.libraries <- function(){
  library(nimble)
  library(MCMCvis)
  library(ggplot2)
  library(bayesplot)
}

load.libraries()

## Load model output ----
#################################

load("model_output/model_output_additive.RData")  # loads: run_js_additive

str(run_js_additive)

## Convergence and distribution check ----
###############################################

# Check traceplots and Rhat for a subset of MCMC iterations:
# Adjust the burn-in cutoff based on inspection of traceplots
burnin_cutoff <- 5001  # example: discard first 5000 iterations
n_total       <- dim(run_js_additive[[1]])[1]

run_additive_list_cutoff <- lapply(run_js_additive, function(x) x[burnin_cutoff:n_total, ])

## Traceplots ----
###############################################

# Intercepts for survival
MCMCtrace(run_additive_list_cutoff,
          params  = c("mean.phi_JD", "mean.phi_JW", "mean.phi_AD", "mean.phi_AW"),
          pdf     = FALSE, ind = TRUE, Rhat = TRUE, n.eff = TRUE,
          iter    = n_total - burnin_cutoff + 1)

# Intercepts for detection
MCMCtrace(run_additive_list_cutoff,
          params  = c("mean.p_JD", "mean.p_JW", "mean.p_AD", "mean.p_AW"),
          pdf     = FALSE, ind = TRUE, Rhat = TRUE, n.eff = TRUE,
          iter    = n_total - burnin_cutoff + 1)

# Beta coefficients for survival
MCMCtrace(run_additive_list_cutoff,
          params  = c("beta.JD.rain_prev", "beta.JD.tmax_current", "beta.JD.dens_current"),
          pdf     = FALSE, ind = TRUE, Rhat = TRUE, n.eff = TRUE,
          iter    = n_total - burnin_cutoff + 1)

MCMCtrace(run_additive_list_cutoff,
          params  = c("beta.JW.rain_current", "beta.JW.tmax_prev", "beta.JW.dens_current"),
          pdf     = FALSE, ind = TRUE, Rhat = TRUE, n.eff = TRUE,
          iter    = n_total - burnin_cutoff + 1)

MCMCtrace(run_additive_list_cutoff,
          params  = c("beta.AD.rain_prev", "beta.AD.tmax_current", "beta.AD.dens_current"),
          pdf     = FALSE, ind = TRUE, Rhat = TRUE, n.eff = TRUE,
          iter    = n_total - burnin_cutoff + 1)

MCMCtrace(run_additive_list_cutoff,
          params  = c("beta.AW.rain_current", "beta.AW.tmax_prev", "beta.AW.dens_current"),
          pdf     = FALSE, ind = TRUE, Rhat = TRUE, n.eff = TRUE,
          iter    = n_total - burnin_cutoff + 1)

# Entry probabilities
MCMCtrace(run_additive_list_cutoff,
          params  = c("gamma"),
          pdf     = FALSE, ind = TRUE, Rhat = TRUE, n.eff = TRUE,
          iter    = n_total - burnin_cutoff + 1)

# Random year effects
MCMCtrace(run_additive_list_cutoff,
          params  = c("eps_phi_JD"),
          pdf     = FALSE, ind = TRUE, Rhat = TRUE, n.eff = TRUE,
          iter    = n_total - burnin_cutoff + 1)

MCMCtrace(run_additive_list_cutoff,
          params  = c("eps_p_JD", "eps_p_JW", "eps_p_AD", "eps_p_AW"),
          pdf     = FALSE, ind = TRUE, Rhat = TRUE, n.eff = TRUE,
          iter    = n_total - burnin_cutoff + 1)


## MCMC summary ----
###############################################

(run_additive_sum_all <- MCMCsummary(run_additive_list_cutoff,
  params = c(
    "mean.phi_JD", "mean.phi_JW", "mean.phi_AD", "mean.phi_AW",
    "mean.p_JD",   "mean.p_JW",   "mean.p_AD",   "mean.p_AW",
    "beta.JD.rain_prev",    "beta.JD.tmax_current", "beta.JD.dens_current",
    "beta.JW.rain_current", "beta.JW.tmax_prev",    "beta.JW.dens_current",
    "beta.AD.rain_prev",    "beta.AD.tmax_current", "beta.AD.dens_current",
    "beta.AW.rain_current", "beta.AW.tmax_prev",    "beta.AW.dens_current",
    "gamma", "pi",
    "sigma_phi_JD", "sigma_phi_JW", "sigma_phi_AD", "sigma_phi_AW",
    "sigma_p_JD",   "sigma_p_JW",   "sigma_p_AD",   "sigma_p_AW",
    "eps_p_JD",     "eps_p_JW",     "eps_p_AD",      "eps_p_AW"
  ),
  round = 3))


## Rhat and effective sample size diagnostics ----
###############################################

# Percentage of Rhat values > 1.1 (convergence threshold):
(pct_rhat_over <- (length(which(run_additive_sum_all$Rhat > 1.1)) /
                   length(run_additive_sum_all$mean)) * 100)

# Distribution of Rhat values:
par(mfrow = c(1, 2))
hist(run_additive_sum_all$Rhat, main = "Rhat distribution", xlab = "Rhat")
abline(v = 1.1, col = "red", lty = 2)

# Effective sample sizes:
plot(run_additive_sum_all$n.eff, ylab = "n.eff", xlab = "Parameter index",
     main = "Effective sample size")
abline(h = 100, col = "red", lty = 2)


## Posterior distribution plots ----
###############################################

# Posterior density plots for beta coefficients (example: JD survival)
# Zero overlap indicates parameter importance

panel_background <- panel_bg(fill = 'white')
color_scheme_set("red")

# JD survival (dry season)
JD_pardist <- mcmc_areas(run_additive_list_cutoff,
                         pars           = c("mean.phi_JD",
                                            "beta.JD.rain_prev",
                                            "beta.JD.tmax_current",
                                            "beta.JD.dens_current"),
                         area_method    = "equal height",
                         prob           = 0.95) +
  panel_background +
  scale_y_discrete(labels = rev(c("Intercept", "Rainfall (prev. wet)",
                                  "Max. Temperature (dry)", "Density"))) +
  theme_classic() +
  theme(axis.text  = element_text(size = 14),
        axis.title = element_text(size = 16))
JD_pardist

# JW survival (wet season)
JW_pardist <- mcmc_areas(run_additive_list_cutoff,
                         pars           = c("mean.phi_JW",
                                            "beta.JW.rain_current",
                                            "beta.JW.tmax_prev",
                                            "beta.JW.dens_current"),
                         area_method    = "equal height",
                         prob           = 0.95) +
  panel_background +
  scale_y_discrete(labels = rev(c("Intercept", "Rainfall (wet)",
                                  "Max. Temperature (prev. dry)", "Density"))) +
  theme_classic() +
  theme(axis.text  = element_text(size = 14),
        axis.title = element_text(size = 16))
JW_pardist

# AD survival (dry season)
AD_pardist <- mcmc_areas(run_additive_list_cutoff,
                         pars           = c("mean.phi_AD",
                                            "beta.AD.rain_prev",
                                            "beta.AD.tmax_current",
                                            "beta.AD.dens_current"),
                         area_method    = "equal height",
                         prob           = 0.95) +
  panel_background +
  scale_y_discrete(labels = rev(c("Intercept", "Rainfall (prev. wet)",
                                  "Max. Temperature (dry)", "Density"))) +
  theme_classic() +
  theme(axis.text  = element_text(size = 14),
        axis.title = element_text(size = 16))
AD_pardist

# AW survival (wet season)
AW_pardist <- mcmc_areas(run_additive_list_cutoff,
                         pars           = c("mean.phi_AW",
                                            "beta.AW.rain_current",
                                            "beta.AW.tmax_prev",
                                            "beta.AW.dens_current"),
                         area_method    = "equal height",
                         prob           = 0.95) +
  panel_background +
  scale_y_discrete(labels = rev(c("Intercept", "Rainfall (wet)",
                                  "Max. Temperature (prev. dry)", "Density"))) +
  theme_classic() +
  theme(axis.text  = element_text(size = 14),
        axis.title = element_text(size = 16))
AW_pardist
