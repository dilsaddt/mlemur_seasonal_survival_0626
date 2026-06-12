############################################################################
# Predictions for models for JAE review

# #######################
# ### ADDITIVE MODELS ###
# #######################

# #########################
# ### FEMALES AND MALES ###
# #########################

# ############################
# ### SEASONAL PREDICTIONS ###
# ############################
# meaning we calculate seasonal predictions like this
# seasonal survival (dry season lasts 5 months, wet season 7 months)
# phi_JD[t] <- phi_month_JD[t]^5
# phi_AD[t] <- phi_month_AD[t]^5
# phi_JW[t] <- phi_month_JW[t]^7
# phi_AW[t] <- phi_month_AW[t]^7

# Not including random year effect in plotting calculations
# so we did not put eps.phi parameters anymore

# I will calculate both 95% and 90% CRIs for all
# then save csv files for each parameter of survival
# also csv files for each recapture
# proper plotting should be in another script

# Date: June 2026
# Authors: Dilsad Dagtekin and Dominik Behr
###########################################################################
## 1. House keeping ----
###########################################################################

rm(list = ls(all = TRUE))

setwd("/path_to_files/")

## 1.1. Load libraries ----
###########################

load.libraries <- function(){
  library(nimble)
  library(parallel)
  library(MCMCvis)
  library(ggplot2)
  library(cowplot)
}

load.libraries()

## 1.2. Load example data ----
#############################################
# Run example_data_for_github.R once to create this file:
load("data/mlemur_example_data_062026.RData")
# Loaded: rain_example, tmax_example, dens_example, yearCov_example,
#         rain_mean, rain_sd, tmax_mean, tmax_sd
seasonal_density <- read.csv("data/mean_seasonal_density.csv")

## 1.3. Prep prediction covariate vectors ----
#############################################

# Smooth prediction vectors over biological ranges (standardized):
rain.stand      <- seq(from = (0    - rain_mean) / rain_sd, to = (1500 - rain_mean) / rain_sd, by = 0.1)
rain.prev.stand <- rain.stand
tmax.stand      <- seq(from = (30   - tmax_mean) / tmax_sd, to = (35   - tmax_mean) / tmax_sd, by = 0.1)
tmax.prev.stand <- tmax.stand

# Density: reconstruct unstandardized series from density model output
df.dens <- seasonal_density
dens    <- df.dens$n.inds.mean
dens    <- c(mean(dens[c(2, 4, 6, 8, 10)]), dens)  # prepend mean of first 5 wet seasons for 1993_2
# Smooth prediction vector over 40–200 individuals:
dens.stand <- seq(from = (40 - mean(dens)) / sd(dens), to = (200 - mean(dens)) / sd(dens), by = 0.1)

# Year covariate (already processed in example data):
yearCov <- yearCov_example


###########################################################################
## 2. Get the model output and prepare the empty slots for results ----
###########################################################################

# Load the most recent additive model output:
model_files <- list.files("model_output/", pattern = "js_model_additive.*\\.RData$",
                           full.names = TRUE)
if(length(model_files) == 0) stop("No additive model output found. Run mlemur_additive_github_0626.R first.")
load(model_files[length(model_files)])
# Combine chains into posterior matrix:
df.post <- do.call(rbind, lapply(run_js_additive, as.matrix))

# Note: example data uses females only.
# To add male predictions: run the model with male data, load that output,
# combine chains into df.post_males, and add to out_list below.
out_list <- list(F = df.post)

# save output csvs here
outdir <- "prediction_output/"

# prep empty lists for results
surv_prob_JD_rain_prev_wet <- list()
surv_prob_JD_rain_prev_wet_add <- list()
surv_prob_JD_tmax_curr_dry <- list()
surv_prob_JD_tmax_curr_dry_add <- list()
surv_prob_JD_dens_curr_dry <- list()
surv_prob_JD_dens_curr_dry_add <- list()


surv_prob_AD_rain_prev_wet <- list()
surv_prob_AD_rain_prev_wet_add <- list()
surv_prob_AD_tmax_curr_dry <-list()
surv_prob_AD_tmax_curr_dry_add <- list()
surv_prob_AD_dens_curr_dry <- list()
surv_prob_AD_dens_curr_dry_add <- list()


surv_prob_JW_rain_current_wet <- list()
surv_prob_JW_rain_current_wet_add <- list()
surv_prob_JW_tmax_prev_dry <- list()
surv_prob_JW_tmax_prev_dry_add <- list()
surv_prob_JW_dens_curr_wet <- list()
surv_prob_JW_dens_curr_wet_add <- list()


surv_prob_AW_rain_current_wet <- list()
surv_prob_AW_rain_current_wet_add <- list()
surv_prob_AW_tmax_prev_dry <- list()
surv_prob_AW_tmax_prev_dry_add <- list()
surv_prob_AW_dens_curr_wet <- list()
surv_prob_AW_dens_curr_wet_add <- list()


# ADDITIVE MODEL ----
###########################################################################
params_add <- colnames(df.post)
params_add

# for both sexes additive model is:
# logit(phi_month_JD[t]) <- mean.phi_JD + beta.JD.rain_prev*rain[t] +
#   beta.JD.tmax_current*tmax[t+1] + 
#   beta.JD.dens_current*dens[t+1] + 
#   beta.JD.tmax_current_dens_current*tmax[t+1]*dens[t+1] +
#   eps_phi_JD[year[t]]

# logit(phi_month_JW[t]) <- mean.phi_JW + beta.JW.rain_current*rain[t+1] +
#   beta.JW.tmax_prev*tmax[t] + 
#   beta.JW.dens_current*dens[t+1] + 
#   beta.JW.tmax_prev_dens_current*tmax[t]*dens[t+1] +
#   eps_phi_JW[year[t]]


# detection: intercept + random year effect

# To get and calculate sex-specific outputs easier

# sex is either M or F; example data uses females only
sex = "F"
out <- out_list[[sex]]
str(out)


# 2.1. DRY SEASON SURVIVALS - JD AND AD ----
###########################################################################

# survival_dry_season: int + rainfall_prev_wet_season + tmax_currenty_dry + density_current

# JD ----------------------------------------------------------------------

mean.phi_JD <- out$mean.phi_JD
beta.JD.rain_prev <- out$beta.JD.rain_prev
beta.JD.tmax_current <- out$beta.JD.tmax_current
beta.JD.dens_current <- out$beta.JD.dens_current

eps_phi_JD <- vector("list", length = length(unique(yearCov$yearCat)))
for (i in 1:length(unique(yearCov$yearCat))) {
  eps_phi_JD[[i]] <- out[,paste0("eps_phi_JD[", i, "]")]
}

###############################################
## RAINFALL OF THE PREVIOUS WET SEASON ----
###############################################

# length.out rain_prev_wet X mcmc samples
surv_JD_rain_prev_wet <- array(NA, dim = c(length(rain.prev.stand), nrow(out)))

for (i in 1:length(rain.prev.stand)){
#  for(t in 1:length(unique(yearCov$yearCat))){
    
    surv_JD_rain_prev_wet[i,] <- plogis(mean.phi_JD 
                                        + beta.JD.rain_prev*rain.prev.stand[i]
                                        + beta.JD.tmax_current*0 # mean tmax value = 0
                                        + beta.JD.dens_current*0)^5 # mean popsize value = 0 and ^5 for seasonal pred of dry season 
                                        # + eps_phi_JD[[t]])^5 # not including random year effect anymore, rest of the script is corrected accordingly
    
 # }
}

str(surv_JD_rain_prev_wet)

# then we take the mean of the mcmc list
pm.surv_JD_rain_prev_wet <- apply(surv_JD_rain_prev_wet, 1, mean)
str(pm.surv_JD_rain_prev_wet)

# then calculate the credible intervals
CRI.surv_JD_rain_prev_wet_95 <- apply(surv_JD_rain_prev_wet, 1, function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_JD_rain_prev_wet_95)
CRI.surv_JD_rain_prev_wet_90 <- apply(surv_JD_rain_prev_wet, 1, function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_JD_rain_prev_wet_90)


surv_prob_JD_rain_prev_wet[[sex]] <- data.frame(state = "JD", 
                                                rain_prev_wet_std = rain.prev.stand,
                                                rain_prev_wet = rain.prev.stand*rain_sd+rain_mean,
                                                param = "survival",
                                                pred = pm.surv_JD_rain_prev_wet,
                                                lower_95 = CRI.surv_JD_rain_prev_wet_95[1,],
                                                upper_95 = CRI.surv_JD_rain_prev_wet_95[2,],
                                                lower_90 = CRI.surv_JD_rain_prev_wet_90[1,],
                                                upper_90 = CRI.surv_JD_rain_prev_wet_90[2,],
                                                sex = sex)


surv_prob_JD_rain_prev_wet[[sex]]

surv_prob_JD_rain_prev_wet_add[[sex]] <- surv_prob_JD_rain_prev_wet[[sex]]

JD_rain_prev_wet_plot <- ggplot(surv_prob_JD_rain_prev_wet_add[[sex]])+
  geom_ribbon(aes(x= rain_prev_wet, ymin = lower_90, ymax = upper_90), fill="blue",alpha=0.3) +
  geom_ribbon(aes(x= rain_prev_wet, ymin = lower_95, ymax = upper_90), fill="lightblue",alpha=0.3) +
  geom_line(aes(rain_prev_wet, pred), col="blue", lwd=3, linetype=1) +
  theme_classic()+
  ylim(0,1) +
  ylab(bquote("Survival Probability " ~ phi ~ "|" ~ .(sex)))+
  xlab('Rainfall of the previous wet season (mm)') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

JD_rain_prev_wet_plot


###############################################
## MAX TEMP OF THE CURRENT DRY SEASON ----
###############################################

# length.out tmax_current_dry X mcmc samples
surv_JD_tmax_curr_dry <- array(NA, dim = c(length(tmax.stand), nrow(out)))

for (i in 1:length(tmax.stand)){
  surv_JD_tmax_curr_dry[i,] <- plogis(mean.phi_JD
                                      + beta.JD.rain_prev*0 # mean rain value = 0
                                      + beta.JD.tmax_current*tmax.stand[i]
                                      + beta.JD.dens_current*0)^5 # mean popsize value = 0 and ^5 for seasonal pred of dry season 
}




str(surv_JD_tmax_curr_dry)

# then we take the mean of the mcmc list
pm.surv_JD_tmax_curr_dry <- apply(surv_JD_tmax_curr_dry, 1, mean)
str(pm.surv_JD_tmax_curr_dry)

# then calculate the credible intervals
CRI.surv_JD_tmax_curr_dry_95 <- apply(surv_JD_tmax_curr_dry, 1, function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_JD_tmax_curr_dry_95)
CRI.surv_JD_tmax_curr_dry_90 <- apply(surv_JD_tmax_curr_dry, 1, function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_JD_tmax_curr_dry_90)


surv_prob_JD_tmax_curr_dry[[sex]] <- data.frame(state = "JD", 
                                                tmax_curr_dry_std = tmax.stand,
                                                tmax_curr_dry = tmax.stand*tmax_sd+tmax_mean,
                                                param = "survival",
                                                pred = pm.surv_JD_tmax_curr_dry,
                                                lower_95 = CRI.surv_JD_tmax_curr_dry_95[1,],
                                                upper_95 = CRI.surv_JD_tmax_curr_dry_95[2,],
                                                lower_90 = CRI.surv_JD_tmax_curr_dry_90[1,],
                                                upper_90 = CRI.surv_JD_tmax_curr_dry_90[2,],
                                                sex = sex)


surv_prob_JD_tmax_curr_dry[[sex]]

surv_prob_JD_tmax_curr_dry_add[[sex]] <- surv_prob_JD_tmax_curr_dry[[sex]]

JD_tmax_curr_dry_plot <- ggplot(surv_prob_JD_tmax_curr_dry_add[[sex]])+
  geom_ribbon(aes(x= tmax_curr_dry, ymin = lower_90, ymax = upper_90), fill="blue",alpha=0.3) +
  geom_ribbon(aes(x= tmax_curr_dry, ymin = lower_95, ymax = upper_90), fill="lightblue",alpha=0.3) +
  geom_line(aes(tmax_curr_dry, pred), col="blue", lwd=3, linetype=1) +
  theme_classic()+
  ylim(0,1) +
  ylab(bquote("Survival Probability " ~ phi ~ "|" ~ .(sex)))+
  xlab('Maximum temperature of the dry season (°C)') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

JD_tmax_curr_dry_plot

###############################################
## POP DENSITY OF THE CURRENT DRY SEASON ----
###############################################

# length.out dens_current X mcmc samples
surv_JD_dens_curr_dry <- array(NA, dim = c(length(dens.stand), nrow(out)))

for (i in 1:length(dens.stand)){
  
  surv_JD_dens_curr_dry[i,] <- plogis(mean.phi_JD 
                                      + beta.JD.rain_prev*0 # mean rain value = 0
                                      + beta.JD.tmax_current*0 # mean tmax value = 0
                                      + beta.JD.dens_current*dens.stand[i])^5 # ^5 for seasonal pred of dry season
  
}




str(surv_JD_dens_curr_dry)

# then we take the mean of the mcmc list
pm.surv_JD_dens_curr_dry <- apply(surv_JD_dens_curr_dry, 1, mean)
str(pm.surv_JD_dens_curr_dry)

# then calculate the credible intervals
CRI.surv_JD_dens_curr_dry_95 <- apply(surv_JD_dens_curr_dry, 1, function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_JD_dens_curr_dry_95)
CRI.surv_JD_dens_curr_dry_90 <- apply(surv_JD_dens_curr_dry, 1, function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_JD_dens_curr_dry_90)


surv_prob_JD_dens_curr_dry[[sex]] <- data.frame(state = "JD", 
                                                dens_curr_dry_std = dens.stand,
                                                dens_curr_dry = dens.stand*sd(dens)+mean(dens),
                                                param = "survival",
                                                pred = pm.surv_JD_dens_curr_dry,
                                                lower_95 = CRI.surv_JD_dens_curr_dry_95[1,],
                                                upper_95 = CRI.surv_JD_dens_curr_dry_95[2,],
                                                lower_90 = CRI.surv_JD_dens_curr_dry_90[1,],
                                                upper_90 = CRI.surv_JD_dens_curr_dry_90[2,],
                                                sex = sex)


surv_prob_JD_dens_curr_dry[[sex]]


surv_prob_JD_dens_curr_dry_add[[sex]] <- surv_prob_JD_dens_curr_dry[[sex]]

JD_dens_curr_dry_plot <- ggplot(surv_prob_JD_dens_curr_dry_add[[sex]])+
  geom_ribbon(aes(x= dens_curr_dry, ymin = lower_90, ymax = upper_90), fill="blue",alpha=0.3) +
  geom_ribbon(aes(x= dens_curr_dry, ymin = lower_95, ymax = upper_90), fill="lightblue",alpha=0.3) +
  geom_line(aes(dens_curr_dry, pred), col="blue", lwd=3, linetype=1) +
  theme_classic()+
  ylim(0,1) +
  ylab(bquote("Survival Probability " ~ phi ~ "|" ~ .(sex)))+
  xlab('Population Density') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

JD_dens_curr_dry_plot


# AD ----------------------------------------------------------------------

mean.phi_AD <- out$mean.phi_AD
beta.AD.rain_prev <- out$beta.AD.rain_prev
beta.AD.tmax_current <- out$beta.AD.tmax_current
beta.AD.dens_current <- out$beta.AD.dens_current

eps_phi_AD <- vector("list", length = length(unique(yearCov$yearCat)))
for (i in 1:length(unique(yearCov$yearCat))) {
  eps_phi_AD[[i]] <- out[,paste0("eps_phi_AD[", i, "]")]
}

###############################################
## RAINFALL OF THE PREVIOUS WET SEASON ----
###############################################

# length.out rain_prev_wet X mcmc samples
surv_AD_rain_prev_wet <- array(NA, dim = c(length(rain.prev.stand), nrow(out)))

for (i in 1:length(rain.prev.stand)){
  #  for(t in 1:length(unique(yearCov$yearCat))){
  
  surv_AD_rain_prev_wet[i,] <- plogis(mean.phi_AD 
                                      + beta.AD.rain_prev*rain.prev.stand[i]
                                      + beta.AD.tmax_current*0 # mean tmax value = 0
                                      + beta.AD.dens_current*0)^5 # mean popsize value = 0 and ^5 for seasonal pred of dry season 
  # + eps_phi_AD[[t]])^5 # not including random year effect anymore, rest of the script is corrected accordingly
  
  # }
}

str(surv_AD_rain_prev_wet)

# then we take the mean of the mcmc list
pm.surv_AD_rain_prev_wet <- apply(surv_AD_rain_prev_wet, 1, mean)
str(pm.surv_AD_rain_prev_wet)

# then calculate the credible intervals
CRI.surv_AD_rain_prev_wet_95 <- apply(surv_AD_rain_prev_wet, 1, function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_AD_rain_prev_wet_95)
CRI.surv_AD_rain_prev_wet_90 <- apply(surv_AD_rain_prev_wet, 1, function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_AD_rain_prev_wet_90)


surv_prob_AD_rain_prev_wet[[sex]] <- data.frame(state = "AD", 
                                                rain_prev_wet_std = rain.prev.stand,
                                                rain_prev_wet = rain.prev.stand*rain_sd+rain_mean,
                                                param = "survival",
                                                pred = pm.surv_AD_rain_prev_wet,
                                                lower_95 = CRI.surv_AD_rain_prev_wet_95[1,],
                                                upper_95 = CRI.surv_AD_rain_prev_wet_95[2,],
                                                lower_90 = CRI.surv_AD_rain_prev_wet_90[1,],
                                                upper_90 = CRI.surv_AD_rain_prev_wet_90[2,],
                                                sex = sex)


surv_prob_AD_rain_prev_wet[[sex]]

surv_prob_AD_rain_prev_wet_add[[sex]] <- surv_prob_AD_rain_prev_wet[[sex]]

AD_rain_prev_wet_plot <- ggplot(surv_prob_AD_rain_prev_wet_add[[sex]])+
  geom_ribbon(aes(x= rain_prev_wet, ymin = lower_90, ymax = upper_90), fill="blue",alpha=0.3) +
  geom_ribbon(aes(x= rain_prev_wet, ymin = lower_95, ymax = upper_90), fill="lightblue",alpha=0.3) +
  geom_line(aes(rain_prev_wet, pred), col="blue", lwd=3, linetype=1) +
  theme_classic()+
  ylim(0,1) +
  ylab(bquote("Survival Probability " ~ phi ~ "|" ~ .(sex)))+
  xlab('Rainfall of the previous wet season (mm)') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

AD_rain_prev_wet_plot


###############################################
## MAX TEMP OF THE CURRENT DRY SEASON ----
###############################################

# length.out tmax_current_dry X mcmc samples
surv_AD_tmax_curr_dry <- array(NA, dim = c(length(tmax.stand), nrow(out)))

for (i in 1:length(tmax.stand)){
  surv_AD_tmax_curr_dry[i,] <- plogis(mean.phi_AD
                                      + beta.AD.rain_prev*0 # mean rain value = 0
                                      + beta.AD.tmax_current*tmax.stand[i]
                                      + beta.AD.dens_current*0)^5 # mean popsize value = 0 and ^5 for seasonal pred of dry season 
}




str(surv_AD_tmax_curr_dry)

# then we take the mean of the mcmc list
pm.surv_AD_tmax_curr_dry <- apply(surv_AD_tmax_curr_dry, 1, mean)
str(pm.surv_AD_tmax_curr_dry)

# then calculate the credible intervals
CRI.surv_AD_tmax_curr_dry_95 <- apply(surv_AD_tmax_curr_dry, 1, function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_AD_tmax_curr_dry_95)
CRI.surv_AD_tmax_curr_dry_90 <- apply(surv_AD_tmax_curr_dry, 1, function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_AD_tmax_curr_dry_90)


surv_prob_AD_tmax_curr_dry[[sex]] <- data.frame(state = "AD", 
                                                tmax_curr_dry_std = tmax.stand,
                                                tmax_curr_dry = tmax.stand*tmax_sd+tmax_mean,
                                                param = "survival",
                                                pred = pm.surv_AD_tmax_curr_dry,
                                                lower_95 = CRI.surv_AD_tmax_curr_dry_95[1,],
                                                upper_95 = CRI.surv_AD_tmax_curr_dry_95[2,],
                                                lower_90 = CRI.surv_AD_tmax_curr_dry_90[1,],
                                                upper_90 = CRI.surv_AD_tmax_curr_dry_90[2,],
                                                sex = sex)


surv_prob_AD_tmax_curr_dry[[sex]]

surv_prob_AD_tmax_curr_dry_add[[sex]] <- surv_prob_AD_tmax_curr_dry[[sex]]

AD_tmax_curr_dry_plot <- ggplot(surv_prob_AD_tmax_curr_dry_add[[sex]])+
  geom_ribbon(aes(x= tmax_curr_dry, ymin = lower_90, ymax = upper_90), fill="blue",alpha=0.3) +
  geom_ribbon(aes(x= tmax_curr_dry, ymin = lower_95, ymax = upper_90), fill="lightblue",alpha=0.3) +
  geom_line(aes(tmax_curr_dry, pred), col="blue", lwd=3, linetype=1) +
  theme_classic()+
  ylim(0,1) +
  ylab(bquote("Survival Probability " ~ phi ~ "|" ~ .(sex)))+
  xlab('Maximum temperature of the dry season (°C)') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

AD_tmax_curr_dry_plot

###############################################
## POP DENSITY OF THE CURRENT DRY SEASON ----
###############################################

# length.out dens_current X mcmc samples
surv_AD_dens_curr_dry <- array(NA, dim = c(length(dens.stand), nrow(out)))

for (i in 1:length(dens.stand)){
  
  surv_AD_dens_curr_dry[i,] <- plogis(mean.phi_AD 
                                      + beta.AD.rain_prev*0 # mean rain value = 0
                                      + beta.AD.tmax_current*0 # mean tmax value = 0
                                      + beta.AD.dens_current*dens.stand[i])^5 # ^5 for seasonal pred of dry season
  
}




str(surv_AD_dens_curr_dry)

# then we take the mean of the mcmc list
pm.surv_AD_dens_curr_dry <- apply(surv_AD_dens_curr_dry, 1, mean)
str(pm.surv_AD_dens_curr_dry)

# then calculate the credible intervals
CRI.surv_AD_dens_curr_dry_95 <- apply(surv_AD_dens_curr_dry, 1, function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_AD_dens_curr_dry_95)
CRI.surv_AD_dens_curr_dry_90 <- apply(surv_AD_dens_curr_dry, 1, function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_AD_dens_curr_dry_90)


surv_prob_AD_dens_curr_dry[[sex]] <- data.frame(state = "AD", 
                                                dens_curr_dry_std = dens.stand,
                                                dens_curr_dry = dens.stand*sd(dens)+mean(dens),
                                                param = "survival",
                                                pred = pm.surv_AD_dens_curr_dry,
                                                lower_95 = CRI.surv_AD_dens_curr_dry_95[1,],
                                                upper_95 = CRI.surv_AD_dens_curr_dry_95[2,],
                                                lower_90 = CRI.surv_AD_dens_curr_dry_90[1,],
                                                upper_90 = CRI.surv_AD_dens_curr_dry_90[2,],
                                                sex = sex)


surv_prob_AD_dens_curr_dry[[sex]]


surv_prob_AD_dens_curr_dry_add[[sex]] <- surv_prob_AD_dens_curr_dry[[sex]]

AD_dens_curr_dry_plot <- ggplot(surv_prob_AD_dens_curr_dry_add[[sex]])+
  geom_ribbon(aes(x= dens_curr_dry, ymin = lower_90, ymax = upper_90), fill="blue",alpha=0.3) +
  geom_ribbon(aes(x= dens_curr_dry, ymin = lower_95, ymax = upper_90), fill="lightblue",alpha=0.3) +
  geom_line(aes(dens_curr_dry, pred), col="blue", lwd=3, linetype=1) +
  theme_classic()+
  ylim(0,1) +
  ylab(bquote("Survival Probability " ~ phi ~ "|" ~ .(sex)))+
  xlab('Population Density') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

AD_dens_curr_dry_plot


# 2.2. WET SEASON SURVIVALS - JW AND AW ----
###########################################################################

# survival_wet_season: int + rainfall_current_wet_season + tmax_prev_season + density_current

# JW ----------------------------------------------------------------------

mean.phi_JW <- out$mean.phi_JW
beta.JW.rain_current <- out$beta.JW.rain_current
beta.JW.tmax_prev <- out$beta.JW.tmax_prev
beta.JW.dens_current <- out$beta.JW.dens_current

eps_phi_JW <- vector("list", length = length(unique(yearCov$yearCat)))
for (i in 1:length(unique(yearCov$yearCat))) {
  eps_phi_JW[[i]] <- out[,paste0("eps_phi_JW[", i, "]")]
}

###############################################
## RAINFALL OF THE CURRENT WET SEASON ----
###############################################

# length.out rain_current_wet X mcmc samples
surv_JW_rain_current_wet <- array(NA, dim = c(length(rain.stand), nrow(out)))

for (i in 1:length(rain.stand)){

    surv_JW_rain_current_wet[i,] <- plogis(mean.phi_JW 
                                           + beta.JW.rain_current*rain.stand[i]
                                           + beta.JW.tmax_prev*0 # mean tmax value = 0
                                           + beta.JW.dens_current*0)^7 # mean popsize value = 0 and ^7 for seasonal pred of wet season
}



str(surv_JW_rain_current_wet)

# then we take the mean of the mcmc list
pm.surv_JW_rain_current_wet <- apply(surv_JW_rain_current_wet, 1, mean)
str(pm.surv_JW_rain_current_wet)

# then calculate the credible intervals
CRI.surv_JW_rain_current_wet_95 <- apply(surv_JW_rain_current_wet, 1, function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_JW_rain_current_wet_95)
CRI.surv_JW_rain_current_wet_90 <- apply(surv_JW_rain_current_wet, 1, function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_JW_rain_current_wet_90)


surv_prob_JW_rain_current_wet[[sex]] <- data.frame(state = "JW", 
                                                   rain_current_wet_std = rain.stand,
                                                   rain_current_wet = rain.stand*rain_sd+rain_mean,
                                                   param = "survival",
                                                   pred = pm.surv_JW_rain_current_wet,
                                                   lower_95 = CRI.surv_JW_rain_current_wet_95[1,],
                                                   upper_95 = CRI.surv_JW_rain_current_wet_95[2,],
                                                   lower_90 = CRI.surv_JW_rain_current_wet_90[1,],
                                                   upper_90 = CRI.surv_JW_rain_current_wet_90[2,],
                                                   sex = sex)


surv_prob_JW_rain_current_wet[[sex]]

surv_prob_JW_rain_current_wet_add[[sex]] <- surv_prob_JW_rain_current_wet[[sex]]

JW_rain_current_wet_plot <- ggplot(surv_prob_JW_rain_current_wet_add[[sex]])+
  geom_ribbon(aes(x= rain_current_wet, ymin = lower_90, ymax = upper_90), fill="blue",alpha=0.3) +
  geom_ribbon(aes(x= rain_current_wet, ymin = lower_95, ymax = upper_90), fill="lightblue",alpha=0.3) +
  geom_line(aes(rain_current_wet, pred), col="blue", lwd=3, linetype=1) +
  theme_classic()+
  ylim(0,1) +
  ylab(bquote("Survival Probability " ~ phi ~ "|" ~ .(sex)))+
  xlab('Rainfall of the current wet season (mm)') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

JW_rain_current_wet_plot


###############################################
## MAX TEMP OF THE PREVIOUS DRY SEASON ----
###############################################

# length.out tmax_prev_dry X mcmc samples
surv_JW_tmax_prev_dry <- array(NA, dim = c(length(tmax.prev.stand), nrow(out)))

for (i in 1:length(tmax.prev.stand)){
  
  surv_JW_tmax_prev_dry[i,] <- plogis(mean.phi_JW 
                                      + beta.JW.rain_current*0 # mean rain value = 0
                                      + beta.JW.tmax_prev*tmax.prev.stand[i]
                                      + beta.JW.dens_current*0)^7 # mean popsize value = 0 and ^7 for seasonal pred of wet season
  
}


str(surv_JW_tmax_prev_dry)

# then we take the mean of the mcmc list
pm.surv_JW_tmax_prev_dry <- apply(surv_JW_tmax_prev_dry, 1, mean)
str(pm.surv_JW_tmax_prev_dry)

# then calculate the credible intervals
CRI.surv_JW_tmax_prev_dry_95 <- apply(surv_JW_tmax_prev_dry, 1, function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_JW_tmax_prev_dry_95)
CRI.surv_JW_tmax_prev_dry_90 <- apply(surv_JW_tmax_prev_dry, 1, function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_JW_tmax_prev_dry_90)


surv_prob_JW_tmax_prev_dry[[sex]] <- data.frame(state = "JW", 
                                                tmax_prev_dry_std = tmax.prev.stand,
                                                tmax_prev_dry = tmax.prev.stand*tmax_sd+tmax_mean,
                                                param = "survival",
                                                pred = pm.surv_JW_tmax_prev_dry,
                                                lower_95 = CRI.surv_JW_tmax_prev_dry_95[1,],
                                                upper_95 = CRI.surv_JW_tmax_prev_dry_95[2,],
                                                lower_90 = CRI.surv_JW_tmax_prev_dry_90[1,],
                                                upper_90 = CRI.surv_JW_tmax_prev_dry_90[2,],
                                                sex = sex)


surv_prob_JW_tmax_prev_dry[[sex]]

surv_prob_JW_tmax_prev_dry_add[[sex]] <- surv_prob_JW_tmax_prev_dry[[sex]]

JW_tmax_prev_dry_plot <- ggplot(surv_prob_JW_tmax_prev_dry_add[[sex]])+
  geom_ribbon(aes(x= tmax_prev_dry, ymin = lower_90, ymax = upper_90), fill="blue",alpha=0.3) +
  geom_ribbon(aes(x= tmax_prev_dry, ymin = lower_95, ymax = upper_90), fill="lightblue",alpha=0.3) +
  geom_line(aes(tmax_prev_dry, pred), col="blue", lwd=3, linetype=1) +
  theme_classic()+
  ylim(0,1) +
  ylab(bquote("Survival Probability " ~ phi ~ "|" ~ .(sex)))+
  xlab('Maximum temperature of the previous dry season (°C)') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

JW_tmax_prev_dry_plot

###############################################
## POP DENSITY OF THE CURRENT WET SEASON ----
###############################################

# length.out dens_current X mcmc samples
surv_JW_dens_curr_wet <- array(NA, dim = c(length(dens.stand), nrow(out)))

for (i in 1:length(dens.stand)){
  
  surv_JW_dens_curr_wet[i,] <- plogis(mean.phi_JW 
                                      + beta.JW.rain_current*0 # mean rain value = 0
                                      + beta.JW.tmax_prev*0 # mean tmax value = 0
                                      + beta.JW.dens_current*dens.stand[i])^7 # ^7 for seasonal pred of wet season
  
}

str(surv_JW_dens_curr_wet)

# then we take the mean of the mcmc list
pm.surv_JW_dens_curr_wet <- apply(surv_JW_dens_curr_wet, 1, mean)
str(pm.surv_JW_dens_curr_wet)

# then calculate the credible intervals
CRI.surv_JW_dens_curr_wet_95 <- apply(surv_JW_dens_curr_wet, 1, function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_JW_dens_curr_wet_95)
CRI.surv_JW_dens_curr_wet_90 <- apply(surv_JW_dens_curr_wet, 1, function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_JW_dens_curr_wet_90)


surv_prob_JW_dens_curr_wet[[sex]] <- data.frame(state = "JW", 
                                                dens_curr_wet_std = dens.stand,
                                                dens_curr_wet = dens.stand*sd(dens)+mean(dens),
                                                param = "survival",
                                                pred = pm.surv_JW_dens_curr_wet,
                                                lower_95 = CRI.surv_JW_dens_curr_wet_95[1,],
                                                upper_95 = CRI.surv_JW_dens_curr_wet_95[2,],
                                                lower_90 = CRI.surv_JW_dens_curr_wet_90[1,],
                                                upper_90 = CRI.surv_JW_dens_curr_wet_90[2,],
                                                sex = sex)


surv_prob_JW_dens_curr_wet[[sex]]


surv_prob_JW_dens_curr_wet_add[[sex]] <- surv_prob_JW_dens_curr_wet[[sex]]

JW_dens_curr_wet_plot <- ggplot(surv_prob_JW_dens_curr_wet_add[[sex]])+
  geom_ribbon(aes(x= dens_curr_wet, ymin = lower_90, ymax = upper_90), fill="blue",alpha=0.3) +
  geom_ribbon(aes(x= dens_curr_wet, ymin = lower_95, ymax = upper_90), fill="lightblue",alpha=0.3) +
  geom_line(aes(dens_curr_wet, pred), col="blue", lwd=3, linetype=1) +
  theme_classic()+
  ylim(0,1) +
  ylab(bquote("Survival Probability " ~ phi ~ "|" ~ .(sex)))+
  xlab('Population Density') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

JW_dens_curr_wet_plot


# AW ----------------------------------------------------------------------

mean.phi_AW <- out$mean.phi_AW
beta.AW.rain_current <- out$beta.AW.rain_current
beta.AW.tmax_prev <- out$beta.AW.tmax_prev
beta.AW.dens_current <- out$beta.AW.dens_current

eps_phi_AW <- vector("list", length = length(unique(yearCov$yearCat)))
for (i in 1:length(unique(yearCov$yearCat))) {
  eps_phi_AW[[i]] <- out[,paste0("eps_phi_AW[", i, "]")]
}

###############################################
## RAINFALL OF THE CURRENT WET SEASON ----
###############################################

# length.out rain_current_wet X mcmc samples
surv_AW_rain_current_wet <- array(NA, dim = c(length(rain.stand), nrow(out)))

for (i in 1:length(rain.stand)){
  
  surv_AW_rain_current_wet[i,] <- plogis(mean.phi_AW 
                                         + beta.AW.rain_current*rain.stand[i]
                                         + beta.AW.tmax_prev*0 # mean tmax value = 0
                                         + beta.AW.dens_current*0)^7 # mean popsize value = 0 and ^7 for seasonal pred of wet season
}



str(surv_AW_rain_current_wet)

# then we take the mean of the mcmc list
pm.surv_AW_rain_current_wet <- apply(surv_AW_rain_current_wet, 1, mean)
str(pm.surv_AW_rain_current_wet)

# then calculate the credible intervals
CRI.surv_AW_rain_current_wet_95 <- apply(surv_AW_rain_current_wet, 1, function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_AW_rain_current_wet_95)
CRI.surv_AW_rain_current_wet_90 <- apply(surv_AW_rain_current_wet, 1, function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_AW_rain_current_wet_90)


surv_prob_AW_rain_current_wet[[sex]] <- data.frame(state = "AW", 
                                                   rain_current_wet_std = rain.stand,
                                                   rain_current_wet = rain.stand*rain_sd+rain_mean,
                                                   param = "survival",
                                                   pred = pm.surv_AW_rain_current_wet,
                                                   lower_95 = CRI.surv_AW_rain_current_wet_95[1,],
                                                   upper_95 = CRI.surv_AW_rain_current_wet_95[2,],
                                                   lower_90 = CRI.surv_AW_rain_current_wet_90[1,],
                                                   upper_90 = CRI.surv_AW_rain_current_wet_90[2,],
                                                   sex = sex)


surv_prob_AW_rain_current_wet[[sex]]

surv_prob_AW_rain_current_wet_add[[sex]] <- surv_prob_AW_rain_current_wet[[sex]]

AW_rain_current_wet_plot <- ggplot(surv_prob_AW_rain_current_wet_add[[sex]])+
  geom_ribbon(aes(x= rain_current_wet, ymin = lower_90, ymax = upper_90), fill="blue",alpha=0.3) +
  geom_ribbon(aes(x= rain_current_wet, ymin = lower_95, ymax = upper_90), fill="lightblue",alpha=0.3) +
  geom_line(aes(rain_current_wet, pred), col="blue", lwd=3, linetype=1) +
  theme_classic()+
  ylim(0,1) +
  ylab(bquote("Survival Probability " ~ phi ~ "|" ~ .(sex)))+
  xlab('Rainfall of the current wet season (mm)') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

AW_rain_current_wet_plot


###############################################
## MAX TEMP OF THE PREVIOUS DRY SEASON ----
###############################################

# length.out tmax_prev_dry X mcmc samples
surv_AW_tmax_prev_dry <- array(NA, dim = c(length(tmax.prev.stand), nrow(out)))

for (i in 1:length(tmax.prev.stand)){
  
  surv_AW_tmax_prev_dry[i,] <- plogis(mean.phi_AW 
                                      + beta.AW.rain_current*0 # mean rain value = 0
                                      + beta.AW.tmax_prev*tmax.prev.stand[i]
                                      + beta.AW.dens_current*0)^7 # mean popsize value = 0 and ^7 for seasonal pred of wet season
  
}


str(surv_AW_tmax_prev_dry)

# then we take the mean of the mcmc list
pm.surv_AW_tmax_prev_dry <- apply(surv_AW_tmax_prev_dry, 1, mean)
str(pm.surv_AW_tmax_prev_dry)

# then calculate the credible intervals
CRI.surv_AW_tmax_prev_dry_95 <- apply(surv_AW_tmax_prev_dry, 1, function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_AW_tmax_prev_dry_95)
CRI.surv_AW_tmax_prev_dry_90 <- apply(surv_AW_tmax_prev_dry, 1, function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_AW_tmax_prev_dry_90)


surv_prob_AW_tmax_prev_dry[[sex]] <- data.frame(state = "AW", 
                                                tmax_prev_dry_std = tmax.prev.stand,
                                                tmax_prev_dry = tmax.prev.stand*tmax_sd+tmax_mean,
                                                param = "survival",
                                                pred = pm.surv_AW_tmax_prev_dry,
                                                lower_95 = CRI.surv_AW_tmax_prev_dry_95[1,],
                                                upper_95 = CRI.surv_AW_tmax_prev_dry_95[2,],
                                                lower_90 = CRI.surv_AW_tmax_prev_dry_90[1,],
                                                upper_90 = CRI.surv_AW_tmax_prev_dry_90[2,],
                                                sex = sex)


surv_prob_AW_tmax_prev_dry[[sex]]

surv_prob_AW_tmax_prev_dry_add[[sex]] <- surv_prob_AW_tmax_prev_dry[[sex]]

AW_tmax_prev_dry_plot <- ggplot(surv_prob_AW_tmax_prev_dry_add[[sex]])+
  geom_ribbon(aes(x= tmax_prev_dry, ymin = lower_90, ymax = upper_90), fill="blue",alpha=0.3) +
  geom_ribbon(aes(x= tmax_prev_dry, ymin = lower_95, ymax = upper_90), fill="lightblue",alpha=0.3) +
  geom_line(aes(tmax_prev_dry, pred), col="blue", lwd=3, linetype=1) +
  theme_classic()+
  ylim(0,1) +
  ylab(bquote("Survival Probability " ~ phi ~ "|" ~ .(sex)))+
  xlab('Maximum temperature of the previous dry season (°C)') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

AW_tmax_prev_dry_plot

###############################################
## POP DENSITY OF THE CURRENT WET SEASON ----
###############################################

# length.out dens_current X mcmc samples
surv_AW_dens_curr_wet <- array(NA, dim = c(length(dens.stand), nrow(out)))

for (i in 1:length(dens.stand)){
  
  surv_AW_dens_curr_wet[i,] <- plogis(mean.phi_AW 
                                      + beta.AW.rain_current*0 # mean rain value = 0
                                      + beta.AW.tmax_prev*0 # mean tmax value = 0
                                      + beta.AW.dens_current*dens.stand[i])^7 # ^7 for seasonal pred of wet season
  
}

str(surv_AW_dens_curr_wet)

# then we take the mean of the mcmc list
pm.surv_AW_dens_curr_wet <- apply(surv_AW_dens_curr_wet, 1, mean)
str(pm.surv_AW_dens_curr_wet)

# then calculate the credible intervals
CRI.surv_AW_dens_curr_wet_95 <- apply(surv_AW_dens_curr_wet, 1, function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_AW_dens_curr_wet_95)
CRI.surv_AW_dens_curr_wet_90 <- apply(surv_AW_dens_curr_wet, 1, function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_AW_dens_curr_wet_90)


surv_prob_AW_dens_curr_wet[[sex]] <- data.frame(state = "AW", 
                                                dens_curr_wet_std = dens.stand,
                                                dens_curr_wet = dens.stand*sd(dens)+mean(dens),
                                                param = "survival",
                                                pred = pm.surv_AW_dens_curr_wet,
                                                lower_95 = CRI.surv_AW_dens_curr_wet_95[1,],
                                                upper_95 = CRI.surv_AW_dens_curr_wet_95[2,],
                                                lower_90 = CRI.surv_AW_dens_curr_wet_90[1,],
                                                upper_90 = CRI.surv_AW_dens_curr_wet_90[2,],
                                                sex = sex)


surv_prob_AW_dens_curr_wet[[sex]]


surv_prob_AW_dens_curr_wet_add[[sex]] <- surv_prob_AW_dens_curr_wet[[sex]]

AW_dens_curr_wet_plot <- ggplot(surv_prob_AW_dens_curr_wet_add[[sex]])+
  geom_ribbon(aes(x= dens_curr_wet, ymin = lower_90, ymax = upper_90), fill="blue",alpha=0.3) +
  geom_ribbon(aes(x= dens_curr_wet, ymin = lower_95, ymax = upper_90), fill="lightblue",alpha=0.3) +
  geom_line(aes(dens_curr_wet, pred), col="blue", lwd=3, linetype=1) +
  theme_classic()+
  ylim(0,1) +
  ylab(bquote("Survival Probability " ~ phi ~ "|" ~ .(sex)))+
  xlab('Population Density') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

AW_dens_curr_wet_plot


###########################################################################
# 3. Save the outputs ----
###########################################################################

# RUN THIS AT THE END WHEN YOU ARE SURE YOU RAN BOTH FEMALE AND MALE RESULTS

# JD

str(surv_prob_JD_rain_prev_wet_add)
surv_prob_JD_rain_prev_wet_add_df <- do.call(rbind, surv_prob_JD_rain_prev_wet_add)
rownames(surv_prob_JD_rain_prev_wet_add_df) <- NULL
surv_prob_JD_rain_prev_wet_add_df$sex_long <- ifelse(surv_prob_JD_rain_prev_wet_add_df$sex == "F", "Female", "Male")
write.csv(surv_prob_JD_rain_prev_wet_add_df,  file = paste0(outdir, "surv_prob_JD_rain_prev_wet_add_df_seasonal_FM.csv"))

str(surv_prob_JD_tmax_curr_dry_add)
surv_prob_JD_tmax_curr_dry_add_df <- do.call(rbind, surv_prob_JD_tmax_curr_dry_add)
rownames(surv_prob_JD_tmax_curr_dry_add_df) <- NULL
surv_prob_JD_tmax_curr_dry_add_df$sex_long <- ifelse(surv_prob_JD_tmax_curr_dry_add_df$sex == "F", "Female", "Male")
write.csv(surv_prob_JD_tmax_curr_dry_add_df,  file = paste0(outdir, "surv_prob_JD_tmax_curr_dry_add_df_seasonal_FM.csv"))

str(surv_prob_JD_dens_curr_dry_add)
surv_prob_JD_dens_curr_dry_add_df <- do.call(rbind, surv_prob_JD_dens_curr_dry_add)
rownames(surv_prob_JD_dens_curr_dry_add_df) <- NULL
surv_prob_JD_dens_curr_dry_add_df$sex_long <- ifelse(surv_prob_JD_dens_curr_dry_add_df$sex == "F", "Female", "Male")
write.csv(surv_prob_JD_dens_curr_dry_add_df,  file = paste0(outdir, "surv_prob_JD_dens_curr_dry_add_df_seasonal_FM.csv"))

# AD


str(surv_prob_AD_rain_prev_wet_add)
surv_prob_AD_rain_prev_wet_add_df <- do.call(rbind, surv_prob_AD_rain_prev_wet_add)
rownames(surv_prob_AD_rain_prev_wet_add_df) <- NULL
surv_prob_AD_rain_prev_wet_add_df$sex_long <- ifelse(surv_prob_AD_rain_prev_wet_add_df$sex == "F", "Female", "Male")
write.csv(surv_prob_AD_rain_prev_wet_add_df,  file = paste0(outdir, "surv_prob_AD_rain_prev_wet_add_df_seasonal_FM.csv"))

str(surv_prob_AD_tmax_curr_dry_add)
surv_prob_AD_tmax_curr_dry_add_df <- do.call(rbind, surv_prob_AD_tmax_curr_dry_add)
rownames(surv_prob_AD_tmax_curr_dry_add_df) <- NULL
surv_prob_AD_tmax_curr_dry_add_df$sex_long <- ifelse(surv_prob_AD_tmax_curr_dry_add_df$sex == "F", "Female", "Male")
write.csv(surv_prob_AD_tmax_curr_dry_add_df,  file = paste0(outdir, "surv_prob_AD_tmax_curr_dry_add_df_seasonal_FM.csv"))

str(surv_prob_AD_dens_curr_dry_add)
surv_prob_AD_dens_curr_dry_add_df <- do.call(rbind, surv_prob_AD_dens_curr_dry_add)
rownames(surv_prob_AD_dens_curr_dry_add_df) <- NULL
surv_prob_AD_dens_curr_dry_add_df$sex_long <- ifelse(surv_prob_AD_dens_curr_dry_add_df$sex == "F", "Female", "Male")
write.csv(surv_prob_AD_dens_curr_dry_add_df,  file = paste0(outdir, "surv_prob_AD_dens_curr_dry_add_df_seasonal_FM.csv"))

# JW

str(surv_prob_JW_rain_current_wet_add)
surv_prob_JW_rain_current_wet_add_df <- do.call(rbind, surv_prob_JW_rain_current_wet_add)
rownames(surv_prob_JW_rain_current_wet_add_df) <- NULL
surv_prob_JW_rain_current_wet_add_df$sex_long <- ifelse(surv_prob_JW_rain_current_wet_add_df$sex == "F", "Female", "Male")
write.csv(surv_prob_JW_rain_current_wet_add_df,  file = paste0(outdir, "surv_prob_JW_rain_current_wet_add_df_seasonal_FM.csv"))

str(surv_prob_JW_tmax_prev_dry_add)
surv_prob_JW_tmax_prev_dry_add_df <- do.call(rbind, surv_prob_JW_tmax_prev_dry_add)
rownames(surv_prob_JW_tmax_prev_dry_add_df) <- NULL
surv_prob_JW_tmax_prev_dry_add_df$sex_long <- ifelse(surv_prob_JW_tmax_prev_dry_add_df$sex == "F", "Female", "Male")
write.csv(surv_prob_JW_tmax_prev_dry_add_df,  file = paste0(outdir, "surv_prob_JW_tmax_prev_dry_add_df_seasonal_FM.csv"))

str(surv_prob_JW_dens_curr_wet_add)
surv_prob_JW_dens_curr_wet_add_df <- do.call(rbind, surv_prob_JW_dens_curr_wet_add)
rownames(surv_prob_JW_dens_curr_wet_add_df) <- NULL
surv_prob_JW_dens_curr_wet_add_df$sex_long <- ifelse(surv_prob_JW_dens_curr_wet_add_df$sex == "F", "Female", "Male")
write.csv(surv_prob_JW_dens_curr_wet_add_df,  file = paste0(outdir, "surv_prob_JW_dens_curr_wet_add_df_seasonal_FM.csv"))

# AW


str(surv_prob_AW_rain_current_wet_add)
surv_prob_AW_rain_current_wet_add_df <- do.call(rbind, surv_prob_AW_rain_current_wet_add)
rownames(surv_prob_AW_rain_current_wet_add_df) <- NULL
surv_prob_AW_rain_current_wet_add_df$sex_long <- ifelse(surv_prob_AW_rain_current_wet_add_df$sex == "F", "Female", "Male")
write.csv(surv_prob_AW_rain_current_wet_add_df,  file = paste0(outdir, "surv_prob_AW_rain_current_wet_add_df_seasonal_FM.csv"))

str(surv_prob_AW_tmax_prev_dry_add)
surv_prob_AW_tmax_prev_dry_add_df <- do.call(rbind, surv_prob_AW_tmax_prev_dry_add)
rownames(surv_prob_AW_tmax_prev_dry_add_df) <- NULL
surv_prob_AW_tmax_prev_dry_add_df$sex_long <- ifelse(surv_prob_AW_tmax_prev_dry_add_df$sex == "F", "Female", "Male")
write.csv(surv_prob_AW_tmax_prev_dry_add_df,  file = paste0(outdir, "surv_prob_AW_tmax_prev_dry_add_df_seasonal_FM.csv"))

str(surv_prob_AW_dens_curr_wet_add)
surv_prob_AW_dens_curr_wet_add_df <- do.call(rbind, surv_prob_AW_dens_curr_wet_add)
rownames(surv_prob_AW_dens_curr_wet_add_df) <- NULL
surv_prob_AW_dens_curr_wet_add_df$sex_long <- ifelse(surv_prob_AW_dens_curr_wet_add_df$sex == "F", "Female", "Male")
write.csv(surv_prob_AW_dens_curr_wet_add_df,  file = paste0(outdir, "surv_prob_AW_dens_curr_wet_add_df_seasonal_FM.csv"))



