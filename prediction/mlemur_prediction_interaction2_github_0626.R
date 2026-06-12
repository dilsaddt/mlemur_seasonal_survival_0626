############################################################################
# Predictions for models for JAE review

# ###########################
# ### INTERACTION MODEL II ###
# ### TMAX:DENS ###
# ###########################
# only getting the interaction effects
# not calculating single effects

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

# I will calculate both 95% and 90% CRIs for all
# then save csv files for each parameter of survival
# also csv files for each recapture
# proper plotting should be in another script

# Not including random year effect in plotting calculations
# so we did not put eps.phi parameters anymore

# Date: June 2026
# Authors: Dilsad Dagtekin and Dominik Behr
###########################################################################

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
df.dens    <- seasonal_density
dens       <- df.dens$n.inds.mean
dens       <- c(mean(dens[c(2, 4, 6, 8, 10)]), dens)  # prepend mean of first 5 wet seasons for 1993_2
dens.stand <- (dens - mean(dens)) / sd(dens)            # standardize actual density values

# for categorical density (high = observed max, low = observed min):
max_dens_std <- dens.stand[which(dens == max(dens))]
min_dens_std <- dens.stand[which(dens == min(dens))]

# Year covariate (already processed in example data):
yearCov <- yearCov_example

###########################################################################
## 2. Get the model output and prepare the empty slots for results ----
###########################################################################

# Load the most recent interaction2 model output:
model_files <- list.files("model_output/", pattern = "js_model_interaction2_tmax_dens.*\\.RData$",
                           full.names = TRUE)
if(length(model_files) == 0) stop("No interaction2 model output found. Run mlemur_interaction2_github_0626.R first.")
load(model_files[length(model_files)])
# Combine chains into posterior matrix:
df.post <- do.call(rbind, lapply(run_js_interaction2, as.matrix))

# Note: example data uses females only.
# To add male predictions: load male model output, combine chains, add to out_list.
out_list <- list(F = df.post)

# save output csvs here
outdir <- "prediction_output/"

# prep empty lists for results
surv_prob_JD_int2_tmax_curr_dry_dens_curr <- list()
surv_prob_JD_int2_tmax_curr_dry_dens_curr_high <- list()
surv_prob_JD_int2_tmax_curr_dry_dens_curr_low <- list()
surv_prob_JD_int2_tmax_curr_dry_dens_curr_cat <- list()

surv_prob_AD_int2_tmax_curr_dry_dens_curr <- list()
surv_prob_AD_int2_tmax_curr_dry_dens_curr_high <- list()
surv_prob_AD_int2_tmax_curr_dry_dens_curr_low <- list()
surv_prob_AD_int2_tmax_curr_dry_dens_curr_cat <- list()

surv_prob_JW_int2_tmax_prev_dry_dens_curr <- list()
surv_prob_JW_int2_tmax_prev_dry_dens_curr_high <- list()
surv_prob_JW_int2_tmax_prev_dry_dens_curr_low <- list()
surv_prob_JW_int2_tmax_prev_dry_dens_curr_cat <- list()

surv_prob_AW_int2_tmax_prev_dry_dens_curr <- list()
surv_prob_AW_int2_tmax_prev_dry_dens_curr_high <- list()
surv_prob_AW_int2_tmax_prev_dry_dens_curr_low <- list()
surv_prob_AW_int2_tmax_prev_dry_dens_curr_cat <- list()

# INTERACTION MODEL II: TMAX:DENS ----
###########################################################################
params_int2 <- colnames(df.post)
params_int2

# for both sexes interaction 2 model is:
# logit(phi_month_JD[t]) <- mean.phi_JD + beta.JD.rain_prev*rain[t] +
#   beta.JD.tmax_current*tmax[t+1] + 
#   beta.JD.dens_current*dens[t+1] + 
#   beta.JD.tmax_current_dens_current*tmax[t+1]*dens[t+1] +
#   eps_phi_JD[year[t]]
# 
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

# survival dry season: int + rainfall_prev_season + tmax_currenty_dry + density_current + random year effect + tmax_currenty_dry:density_current

# JD ----------------------------------------------------------------------

mean.phi_JD <- out$mean.phi_JD
beta.JD.rain_prev <- out$beta.JD.rain_prev
beta.JD.tmax_current <- out$beta.JD.tmax_current
beta.JD.dens_current <- out$beta.JD.dens_current
beta.JD.tmax_current_dens_current <- out$beta.JD.tmax_current_dens_current

eps_phi_JD <- vector("list", length = length(unique(yearCov$yearCat)))
for (i in 1:length(unique(yearCov$yearCat))) {
  eps_phi_JD[[i]] <- out[,paste0("eps_phi_JD[", i, "]")]
}

###################################################################
## TMAX OF THE CURRENT DRY SEASON & POPDENS INTERACTION ----
###################################################################
# for the output data frames below
tmax_curr_dry_std = tmax.stand
tmax_curr_dry = tmax.stand*tmax_sd+tmax_mean

# length.out rain_prev_wet X length.oyt.dens_current X mcmc samples
surv_JD_tmax_curr_dry_dens_curr <- array(NA, dim = c(length(tmax.stand), length(dens.stand), nrow(out)))

for (i in 1:length(tmax.stand)){
  for (j in 1: length(dens.stand)){
    #  for(t in 1:length(unique(yearCov$yearCat))){
    
    surv_JD_tmax_curr_dry_dens_curr[i,j,] <- plogis(mean.phi_JD 
                                                    + beta.JD.rain_prev*0 # mean tmax value = 0
                                                    + beta.JD.tmax_current*tmax.stand[i] # mean tmax value = 0
                                                    + beta.JD.dens_current*dens.stand[j]
                                                    + beta.JD.tmax_current_dens_current*tmax.stand[i]*dens.stand[j])^5 # ^5 for seasonal pred of dry season
    # + eps_phi_JD[[t]])^5 # not including random year effect anymore, rest of the script is corrected accordingly 
    
    # }
  }
}

str(surv_JD_tmax_curr_dry_dens_curr)

# then we take the mean of the mcmc list
pm.surv_JD_tmax_curr_dry_dens_curr <- apply(surv_JD_tmax_curr_dry_dens_curr, c(1,2), mean)
str(pm.surv_JD_tmax_curr_dry_dens_curr)

# then calculate the credible intervals
CRI.surv_JD_tmax_curr_dry_dens_curr_95 <- apply(surv_JD_tmax_curr_dry_dens_curr, c(1,2), function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_JD_tmax_curr_dry_dens_curr_95)
CRI.surv_JD_tmax_curr_dry_dens_curr_90 <- apply(surv_JD_tmax_curr_dry_dens_curr, c(1,2), function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_JD_tmax_curr_dry_dens_curr_90)


surv_prob_JD_int2_tmax_curr_dry_dens_curr[[sex]] <- expand.grid(tmax_std = tmax_curr_dry_std, dens_curr_std = dens.stand)


surv_prob_JD_int2_tmax_curr_dry_dens_curr[[sex]]$tmax_curr = rep(round(tmax_curr_dry, digits = 2), length = length(tmax.stand))
surv_prob_JD_int2_tmax_curr_dry_dens_curr[[sex]]$dens_curr = rep(round(dens), each = length(tmax.stand))

surv_prob_JD_int2_tmax_curr_dry_dens_curr[[sex]]$state = "JD"
surv_prob_JD_int2_tmax_curr_dry_dens_curr[[sex]]$param = "survival"
surv_prob_JD_int2_tmax_curr_dry_dens_curr[[sex]]$pred = c(pm.surv_JD_tmax_curr_dry_dens_curr)
surv_prob_JD_int2_tmax_curr_dry_dens_curr[[sex]]$lower_95 = c(CRI.surv_JD_tmax_curr_dry_dens_curr_95[1,,])
surv_prob_JD_int2_tmax_curr_dry_dens_curr[[sex]]$upper_95 = c(CRI.surv_JD_tmax_curr_dry_dens_curr_95[2,,])
surv_prob_JD_int2_tmax_curr_dry_dens_curr[[sex]]$lower_90 = c(CRI.surv_JD_tmax_curr_dry_dens_curr_90[1,,])
surv_prob_JD_int2_tmax_curr_dry_dens_curr[[sex]]$upper_90 = c(CRI.surv_JD_tmax_curr_dry_dens_curr_90[2,,])
surv_prob_JD_int2_tmax_curr_dry_dens_curr[[sex]]$sex = sex

surv_prob_JD_int2_tmax_curr_dry_dens_curr[[sex]]


JD_tmax_dens_int2_plot <- ggplot(surv_prob_JD_int2_tmax_curr_dry_dens_curr[[sex]])+
  geom_ribbon(aes(x= tmax_curr, ymin = lower_90, ymax = upper_90, group = dens_curr , fill = dens_curr ), alpha=0.03) +
  geom_line(aes(tmax_curr, pred, group = dens_curr, col = dens_curr ), lwd=3, linetype=1) +
  theme_classic()+
  scale_fill_continuous(type = "gradient", low = "#bcbddc", high = "#3f007d") +
  scale_color_continuous(type = "gradient", low = "#bcbddc", high = "#3f007d") +
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

JD_tmax_dens_int2_plot

##############################################################################
## TMAX OF THE CURRENT DRY SEASON & POPDENS INTERACTION CATEGORICAL ----
##############################################################################
# for the output data frames below
tmax_curr_dry_std = tmax.stand
tmax_curr_dry = tmax.stand*tmax_sd+tmax_mean

# FOR HIGH DENSITY = 194 INDS, dens_std = 3.155485
# FOR LOW DENSITY = 40 INDS, dens_std = -1.322384

### HIGH DENSITY ###

# length.out rain_prev_wet X mcmc samples
surv_JD_tmax_curr_dry_dens_curr_high <- array(NA, dim = c(length(tmax.stand), nrow(out)))

for (i in 1:length(tmax.stand)){
  
  surv_JD_tmax_curr_dry_dens_curr_high[i,] <- plogis(mean.phi_JD 
                                                     + beta.JD.rain_prev*0 # mean tmax value = 0
                                                     + beta.JD.tmax_current*tmax.stand[i] # mean tmax value = 0
                                                     + beta.JD.dens_current*max_dens_std
                                                     + beta.JD.tmax_current_dens_current*tmax.stand[i]*max_dens_std)^5 # ^5 for seasonal pred of dry season
  
}

str(surv_JD_tmax_curr_dry_dens_curr_high)

# then we take the mean of the mcmc list
pm.surv_JD_tmax_curr_dry_dens_curr_high <- apply(surv_JD_tmax_curr_dry_dens_curr_high, 1, mean)
str(pm.surv_JD_tmax_curr_dry_dens_curr_high)

# then calculate the credible intervals
CRI.surv_JD_tmax_curr_dry_dens_curr_high_95 <- apply(surv_JD_tmax_curr_dry_dens_curr_high, 1, function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_JD_tmax_curr_dry_dens_curr_high_95)
CRI.surv_JD_tmax_curr_dry_dens_curr_high_90 <- apply(surv_JD_tmax_curr_dry_dens_curr_high, 1, function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_JD_tmax_curr_dry_dens_curr_high_90)


surv_prob_JD_int2_tmax_curr_dry_dens_curr_high[[sex]] <- data.frame(state = "JD", 
                                                                    tmax_curr_std = tmax_curr_dry_std,
                                                                    tmax_curr = round(tmax_curr_dry, digits = 2),
                                                                    dens_curr_std = max_dens_std,
                                                                    dens_curr = round(max(dens)),
                                                                    dens_cat = "High Density",
                                                                    param = "survival",
                                                                    pred = pm.surv_JD_tmax_curr_dry_dens_curr_high,
                                                                    lower_95 = CRI.surv_JD_tmax_curr_dry_dens_curr_high_95[1,],
                                                                    upper_95 = CRI.surv_JD_tmax_curr_dry_dens_curr_high_95[2,],
                                                                    lower_90 = CRI.surv_JD_tmax_curr_dry_dens_curr_high_90[1,],
                                                                    upper_90 = CRI.surv_JD_tmax_curr_dry_dens_curr_high_90[2,],
                                                                    sex = sex)


surv_prob_JD_int2_tmax_curr_dry_dens_curr_high[[sex]]


### LOW DENSITY ###

# length.out rain_prev_wet X mcmc samples
surv_JD_tmax_curr_dry_dens_curr_low <- array(NA, dim = c(length(tmax.stand), nrow(out)))

for (i in 1:length(tmax.stand)){
  
  surv_JD_tmax_curr_dry_dens_curr_low[i,] <- plogis(mean.phi_JD 
                                                    + beta.JD.rain_prev*0 # mean tmax value = 0
                                                    + beta.JD.tmax_current*tmax.stand[i] # mean tmax value = 0
                                                    + beta.JD.dens_current*min_dens_std
                                                    + beta.JD.tmax_current_dens_current*tmax.stand[i]*min_dens_std)^5 # ^5 for seasonal pred of dry season
  
}

str(surv_JD_tmax_curr_dry_dens_curr_low)

# then we take the mean of the mcmc list
pm.surv_JD_tmax_curr_dry_dens_curr_low <- apply(surv_JD_tmax_curr_dry_dens_curr_low, 1, mean)
str(pm.surv_JD_tmax_curr_dry_dens_curr_low)

# then calculate the credible intervals
CRI.surv_JD_tmax_curr_dry_dens_curr_low_95 <- apply(surv_JD_tmax_curr_dry_dens_curr_low, 1, function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_JD_tmax_curr_dry_dens_curr_low_95)
CRI.surv_JD_tmax_curr_dry_dens_curr_low_90 <- apply(surv_JD_tmax_curr_dry_dens_curr_low, 1, function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_JD_tmax_curr_dry_dens_curr_low_90)


surv_prob_JD_int2_tmax_curr_dry_dens_curr_low[[sex]] <- data.frame(state = "JD", 
                                                                   tmax_curr_std = tmax_curr_dry_std,
                                                                   tmax_curr = round(tmax_curr_dry, digits = 2),
                                                                   dens_curr_std = min_dens_std,
                                                                   dens_curr = round(min(dens)),
                                                                   dens_cat = "Low Density",
                                                                   param = "survival",
                                                                   pred = pm.surv_JD_tmax_curr_dry_dens_curr_low,
                                                                   lower_95 = CRI.surv_JD_tmax_curr_dry_dens_curr_low_95[1,],
                                                                   upper_95 = CRI.surv_JD_tmax_curr_dry_dens_curr_low_95[2,],
                                                                   lower_90 = CRI.surv_JD_tmax_curr_dry_dens_curr_low_90[1,],
                                                                   upper_90 = CRI.surv_JD_tmax_curr_dry_dens_curr_low_90[2,],
                                                                   sex = sex)


surv_prob_JD_int2_tmax_curr_dry_dens_curr_low[[sex]]


surv_prob_JD_int2_tmax_curr_dry_dens_curr_cat[[sex]] <- rbind(surv_prob_JD_int2_tmax_curr_dry_dens_curr_low[[sex]], surv_prob_JD_int2_tmax_curr_dry_dens_curr_high[[sex]])

surv_JD_int2_plot <- ggplot(surv_prob_JD_int2_tmax_curr_dry_dens_curr_cat[[sex]])+
  geom_ribbon(aes(x= tmax_curr, ymin = lower_90, ymax = upper_90), fill="blue",alpha=0.3) +
  geom_ribbon(aes(x= tmax_curr, ymin = lower_95, ymax = upper_90), fill="lightblue",alpha=0.3) +
  geom_line(aes(tmax_curr, pred), col="blue", lwd=3, linetype=1) +
  theme_classic()+
  facet_grid(~dens_cat) + 
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

surv_JD_int2_plot


# AD ----------------------------------------------------------------------

mean.phi_AD <- out$mean.phi_AD
beta.AD.rain_prev <- out$beta.AD.rain_prev
beta.AD.tmax_current <- out$beta.AD.tmax_current
beta.AD.dens_current <- out$beta.AD.dens_current
beta.AD.tmax_current_dens_current <- out$beta.AD.tmax_current_dens_current

eps_phi_AD <- vector("list", length = length(unique(yearCov$yearCat)))
for (i in 1:length(unique(yearCov$yearCat))) {
  eps_phi_AD[[i]] <- out[,paste0("eps_phi_AD[", i, "]")]
}

###################################################################
## TMAX OF THE CURRENT DRY SEASON & POPDENS INTERACTION ----
###################################################################
# for the output data frames below
tmax_curr_dry_std = tmax.stand
tmax_curr_dry = tmax.stand*tmax_sd+tmax_mean

# length.out rain_prev_wet X length.oyt.dens_current X mcmc samples
surv_AD_tmax_curr_dry_dens_curr <- array(NA, dim = c(length(tmax.stand), length(dens.stand), nrow(out)))

for (i in 1:length(tmax.stand)){
  for (j in 1: length(dens.stand)){
    
    surv_AD_tmax_curr_dry_dens_curr[i,j,] <- plogis(mean.phi_AD 
                                                    + beta.AD.rain_prev*0 # mean tmax value = 0
                                                    + beta.AD.tmax_current*tmax.stand[i] # mean tmax value = 0
                                                    + beta.AD.dens_current*dens.stand[j]
                                                    + beta.AD.tmax_current_dens_current*tmax.stand[i]*dens.stand[j])^5 # ^5 for seasonal pred of dry season
  }
}

str(surv_AD_tmax_curr_dry_dens_curr)

# then we take the mean of the mcmc list
pm.surv_AD_tmax_curr_dry_dens_curr <- apply(surv_AD_tmax_curr_dry_dens_curr, c(1,2), mean)
str(pm.surv_AD_tmax_curr_dry_dens_curr)

# then calculate the credible intervals
CRI.surv_AD_tmax_curr_dry_dens_curr_95 <- apply(surv_AD_tmax_curr_dry_dens_curr, c(1,2), function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_AD_tmax_curr_dry_dens_curr_95)
CRI.surv_AD_tmax_curr_dry_dens_curr_90 <- apply(surv_AD_tmax_curr_dry_dens_curr, c(1,2), function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_AD_tmax_curr_dry_dens_curr_90)


surv_prob_AD_int2_tmax_curr_dry_dens_curr[[sex]] <- expand.grid(tmax_std = tmax_curr_dry_std, dens_curr_std = dens.stand)


surv_prob_AD_int2_tmax_curr_dry_dens_curr[[sex]]$tmax_curr = rep(round(tmax_curr_dry, digits = 2), length = length(tmax.stand))
surv_prob_AD_int2_tmax_curr_dry_dens_curr[[sex]]$dens_curr = rep(round(dens), each = length(tmax.stand))

surv_prob_AD_int2_tmax_curr_dry_dens_curr[[sex]]$state = "AD"
surv_prob_AD_int2_tmax_curr_dry_dens_curr[[sex]]$param = "survival"
surv_prob_AD_int2_tmax_curr_dry_dens_curr[[sex]]$pred = c(pm.surv_AD_tmax_curr_dry_dens_curr)
surv_prob_AD_int2_tmax_curr_dry_dens_curr[[sex]]$lower_95 = c(CRI.surv_AD_tmax_curr_dry_dens_curr_95[1,,])
surv_prob_AD_int2_tmax_curr_dry_dens_curr[[sex]]$upper_95 = c(CRI.surv_AD_tmax_curr_dry_dens_curr_95[2,,])
surv_prob_AD_int2_tmax_curr_dry_dens_curr[[sex]]$lower_90 = c(CRI.surv_AD_tmax_curr_dry_dens_curr_90[1,,])
surv_prob_AD_int2_tmax_curr_dry_dens_curr[[sex]]$upper_90 = c(CRI.surv_AD_tmax_curr_dry_dens_curr_90[2,,])
surv_prob_AD_int2_tmax_curr_dry_dens_curr[[sex]]$sex = sex

surv_prob_AD_int2_tmax_curr_dry_dens_curr[[sex]]


AD_tmax_dens_int2_plot <- ggplot(surv_prob_AD_int2_tmax_curr_dry_dens_curr[[sex]])+
  geom_ribbon(aes(x= tmax_curr, ymin = lower_90, ymax = upper_90, group = dens_curr , fill = dens_curr ), alpha=0.03) +
  geom_line(aes(tmax_curr, pred, group = dens_curr, col = dens_curr ), lwd=3, linetype=1) +
  theme_classic()+
  scale_fill_continuous(type = "gradient", low = "#bcbddc", high = "#3f007d") +
  scale_color_continuous(type = "gradient", low = "#bcbddc", high = "#3f007d") +
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

AD_tmax_dens_int2_plot

##############################################################################
## TMAX OF THE CURRENT DRY SEASON & POPDENS INTERACTION CATEGORICAL ----
##############################################################################
# for the output data frames below
tmax_curr_dry_std = tmax.stand
tmax_curr_dry = tmax.stand*tmax_sd+tmax_mean

# FOR HIGH DENSITY = 194 INDS, dens_std = 3.155485
# FOR LOW DENSITY = 40 INDS, dens_std = -1.322384

### HIGH DENSITY ###

# length.out rain_prev_wet X mcmc samples
surv_AD_tmax_curr_dry_dens_curr_high <- array(NA, dim = c(length(tmax.stand), nrow(out)))

for (i in 1:length(tmax.stand)){
  
  surv_AD_tmax_curr_dry_dens_curr_high[i,] <- plogis(mean.phi_AD 
                                                     + beta.AD.rain_prev*0 # mean tmax value = 0
                                                     + beta.AD.tmax_current*tmax.stand[i] # mean tmax value = 0
                                                     + beta.AD.dens_current*max_dens_std
                                                     + beta.AD.tmax_current_dens_current*tmax.stand[i]*max_dens_std)^5 # ^5 for seasonal pred of dry season
  
}

str(surv_AD_tmax_curr_dry_dens_curr_high)

# then we take the mean of the mcmc list
pm.surv_AD_tmax_curr_dry_dens_curr_high <- apply(surv_AD_tmax_curr_dry_dens_curr_high, 1, mean)
str(pm.surv_AD_tmax_curr_dry_dens_curr_high)

# then calculate the credible intervals
CRI.surv_AD_tmax_curr_dry_dens_curr_high_95 <- apply(surv_AD_tmax_curr_dry_dens_curr_high, 1, function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_AD_tmax_curr_dry_dens_curr_high_95)
CRI.surv_AD_tmax_curr_dry_dens_curr_high_90 <- apply(surv_AD_tmax_curr_dry_dens_curr_high, 1, function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_AD_tmax_curr_dry_dens_curr_high_90)


surv_prob_AD_int2_tmax_curr_dry_dens_curr_high[[sex]] <- data.frame(state = "AD", 
                                                                    tmax_curr_std = tmax_curr_dry_std,
                                                                    tmax_curr = round(tmax_curr_dry, digits = 2),
                                                                    dens_curr_std = max_dens_std,
                                                                    dens_curr = round(max(dens)),
                                                                    dens_cat = "High Density",
                                                                    param = "survival",
                                                                    pred = pm.surv_AD_tmax_curr_dry_dens_curr_high,
                                                                    lower_95 = CRI.surv_AD_tmax_curr_dry_dens_curr_high_95[1,],
                                                                    upper_95 = CRI.surv_AD_tmax_curr_dry_dens_curr_high_95[2,],
                                                                    lower_90 = CRI.surv_AD_tmax_curr_dry_dens_curr_high_90[1,],
                                                                    upper_90 = CRI.surv_AD_tmax_curr_dry_dens_curr_high_90[2,],
                                                                    sex = sex)


surv_prob_AD_int2_tmax_curr_dry_dens_curr_high[[sex]]


### LOW DENSITY ###

# length.out rain_prev_wet X mcmc samples
surv_AD_tmax_curr_dry_dens_curr_low <- array(NA, dim = c(length(tmax.stand), nrow(out)))

for (i in 1:length(tmax.stand)){
  
  surv_AD_tmax_curr_dry_dens_curr_low[i,] <- plogis(mean.phi_AD 
                                                    + beta.AD.rain_prev*0 # mean tmax value = 0
                                                    + beta.AD.tmax_current*tmax.stand[i] # mean tmax value = 0
                                                    + beta.AD.dens_current*min_dens_std
                                                    + beta.AD.tmax_current_dens_current*tmax.stand[i]*min_dens_std)^5 # ^5 for seasonal pred of dry season
  
}

str(surv_AD_tmax_curr_dry_dens_curr_low)

# then we take the mean of the mcmc list
pm.surv_AD_tmax_curr_dry_dens_curr_low <- apply(surv_AD_tmax_curr_dry_dens_curr_low, 1, mean)
str(pm.surv_AD_tmax_curr_dry_dens_curr_low)

# then calculate the credible intervals
CRI.surv_AD_tmax_curr_dry_dens_curr_low_95 <- apply(surv_AD_tmax_curr_dry_dens_curr_low, 1, function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_AD_tmax_curr_dry_dens_curr_low_95)
CRI.surv_AD_tmax_curr_dry_dens_curr_low_90 <- apply(surv_AD_tmax_curr_dry_dens_curr_low, 1, function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_AD_tmax_curr_dry_dens_curr_low_90)


surv_prob_AD_int2_tmax_curr_dry_dens_curr_low[[sex]] <- data.frame(state = "AD", 
                                                                   tmax_curr_std = tmax_curr_dry_std,
                                                                   tmax_curr = round(tmax_curr_dry, digits = 2),
                                                                   dens_curr_std = min_dens_std,
                                                                   dens_curr = round(min(dens)),
                                                                   dens_cat = "Low Density",
                                                                   param = "survival",
                                                                   pred = pm.surv_AD_tmax_curr_dry_dens_curr_low,
                                                                   lower_95 = CRI.surv_AD_tmax_curr_dry_dens_curr_low_95[1,],
                                                                   upper_95 = CRI.surv_AD_tmax_curr_dry_dens_curr_low_95[2,],
                                                                   lower_90 = CRI.surv_AD_tmax_curr_dry_dens_curr_low_90[1,],
                                                                   upper_90 = CRI.surv_AD_tmax_curr_dry_dens_curr_low_90[2,],
                                                                   sex = sex)


surv_prob_AD_int2_tmax_curr_dry_dens_curr_low[[sex]]


surv_prob_AD_int2_tmax_curr_dry_dens_curr_cat[[sex]] <- rbind(surv_prob_AD_int2_tmax_curr_dry_dens_curr_low[[sex]], surv_prob_AD_int2_tmax_curr_dry_dens_curr_high[[sex]])

surv_AD_int2_plot <- ggplot(surv_prob_AD_int2_tmax_curr_dry_dens_curr_cat[[sex]])+
  geom_ribbon(aes(x= tmax_curr, ymin = lower_90, ymax = upper_90), fill="blue",alpha=0.3) +
  geom_ribbon(aes(x= tmax_curr, ymin = lower_95, ymax = upper_90), fill="lightblue",alpha=0.3) +
  geom_line(aes(tmax_curr, pred), col="blue", lwd=3, linetype=1) +
  theme_classic()+
  facet_grid(~dens_cat) + 
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

surv_AD_int2_plot

# 2.2. WET SEASON SURVIVALS - JW AND AW ----
###########################################################################

# survival wet season: int + rainfall_current_Wet_season + tmax_prev_season + density_current + random year effect + tmax_prev_season:density_current

# JW ----------------------------------------------------------------------
mean.phi_JW <- out$mean.phi_JW
beta.JW.rain_current <- out$beta.JW.rain_current
beta.JW.tmax_prev <- out$beta.JW.tmax_prev
beta.JW.dens_current <- out$beta.JW.dens_current
beta.JW.tmax_prev_dens_current <- out$beta.JW.tmax_prev_dens_current

eps_phi_JW <- vector("list", length = length(unique(yearCov$yearCat)))
for (i in 1:length(unique(yearCov$yearCat))) {
  eps_phi_JW[[i]] <- out[,paste0("eps_phi_JW[", i, "]")]
}

###################################################################
## TMAX OF THE PREVIOUS DRY SEASON & POPDENS INTERACTION ----
###################################################################
# for the output data frames below
tmax_prev_dry_std = tmax.prev.stand
tmax_prev_dry = tmax.prev.stand*tmax_sd+tmax_mean

# length.out rain_prev_wet X length.oyt.dens_current X mcmc samples
surv_JW_tmax_prev_dry_dens_curr <- array(NA, dim = c(length(tmax.prev.stand), length(dens.stand), nrow(out)))

for (i in 1:length(tmax.prev.stand)){
  for (j in 1: length(dens.stand)){
    
    surv_JW_tmax_prev_dry_dens_curr[i,j,] <- plogis(mean.phi_JW 
                                                    + beta.JW.rain_current*0 # mean rain value = 0
                                                    + beta.JW.tmax_prev*tmax.prev.stand[i]
                                                    + beta.JW.dens_current*dens.stand[j]
                                                    + beta.JW.tmax_prev_dens_current*tmax.prev.stand[i]*dens.stand[j])^7 # ^7 for seasonal pred of wet season
  }
}

str(surv_JW_tmax_prev_dry_dens_curr)

# then we take the mean of the mcmc list
pm.surv_JW_tmax_prev_dry_dens_curr <- apply(surv_JW_tmax_prev_dry_dens_curr, c(1,2), mean)
str(pm.surv_JW_tmax_prev_dry_dens_curr)

# then calculate the credible intervals
CRI.surv_JW_tmax_prev_dry_dens_curr_95 <- apply(surv_JW_tmax_prev_dry_dens_curr, c(1,2), function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_JW_tmax_prev_dry_dens_curr_95)
CRI.surv_JW_tmax_prev_dry_dens_curr_90 <- apply(surv_JW_tmax_prev_dry_dens_curr, c(1,2), function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_JW_tmax_prev_dry_dens_curr_90)


surv_prob_JW_int2_tmax_prev_dry_dens_curr[[sex]] <- expand.grid(tmax_prev_dry_std = tmax_prev_dry_std, dens_curr_std = dens.stand)


surv_prob_JW_int2_tmax_prev_dry_dens_curr[[sex]]$tmax_prev_dry = rep(round(tmax_prev_dry, digits = 2), length = length(tmax.prev.stand))
surv_prob_JW_int2_tmax_prev_dry_dens_curr[[sex]]$dens_curr = rep(round(dens), each = length(tmax.prev.stand))

surv_prob_JW_int2_tmax_prev_dry_dens_curr[[sex]]$state = "JW"
surv_prob_JW_int2_tmax_prev_dry_dens_curr[[sex]]$param = "survival"
surv_prob_JW_int2_tmax_prev_dry_dens_curr[[sex]]$pred = c(pm.surv_JW_tmax_prev_dry_dens_curr)
surv_prob_JW_int2_tmax_prev_dry_dens_curr[[sex]]$lower_95 = c(CRI.surv_JW_tmax_prev_dry_dens_curr_95[1,,])
surv_prob_JW_int2_tmax_prev_dry_dens_curr[[sex]]$upper_95 = c(CRI.surv_JW_tmax_prev_dry_dens_curr_95[2,,])
surv_prob_JW_int2_tmax_prev_dry_dens_curr[[sex]]$lower_90 = c(CRI.surv_JW_tmax_prev_dry_dens_curr_90[1,,])
surv_prob_JW_int2_tmax_prev_dry_dens_curr[[sex]]$upper_90 = c(CRI.surv_JW_tmax_prev_dry_dens_curr_90[2,,])
surv_prob_JW_int2_tmax_prev_dry_dens_curr[[sex]]$sex = sex

surv_prob_JW_int2_tmax_prev_dry_dens_curr[[sex]]


JW_tmax_dens_int2_plot <- ggplot(surv_prob_JW_int2_tmax_prev_dry_dens_curr[[sex]])+
  geom_ribbon(aes(x= tmax_prev_dry, ymin = lower_90, ymax = upper_90, group = dens_curr , fill = dens_curr ), alpha=0.03) +
  geom_line(aes(tmax_prev_dry, pred, group = dens_curr, col = dens_curr ), lwd=3, linetype=1) +
  theme_classic()+
  scale_fill_continuous(type = "gradient", low = "#bcbddc", high = "#3f007d") +
  scale_color_continuous(type = "gradient", low = "#bcbddc", high = "#3f007d") +
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

JW_tmax_dens_int2_plot

##############################################################################
## TMAX OF THE PREVIOUS DRY SEASON & POPDENS INTERACTION CATEGORICAL ----
##############################################################################
# FOR HIGH DENSITY = 194 INDS, dens_std = 3.155485
# FOR LOW DENSITY = 40 INDS, dens_std = -1.322384

# for the output data frames below
tmax_prev_dry_std = tmax.prev.stand
tmax_prev_dry = tmax.prev.stand*tmax_sd+tmax_mean

### HIGH DENSITY ###

# length.out rain_prev_wet X mcmc samples
surv_JW_tmax_prev_dry_dens_curr_high <- array(NA, dim = c(length(tmax.prev.stand), nrow(out)))

for (i in 1:length(tmax.prev.stand)){
  
  surv_JW_tmax_prev_dry_dens_curr_high[i,] <- plogis(mean.phi_JW 
                                                     + beta.JW.rain_current*0 # mean rain value = 0
                                                     + beta.JW.tmax_prev*tmax.prev.stand[i]
                                                     + beta.JW.dens_current*max_dens_std
                                                     + beta.JW.tmax_prev_dens_current*tmax.prev.stand[i]*max_dens_std)^7 # ^7 for seasonal pred of wet season
  
}

str(surv_JW_tmax_prev_dry_dens_curr_high)

# then we take the mean of the mcmc list
pm.surv_JW_tmax_prev_dry_dens_curr_high <- apply(surv_JW_tmax_prev_dry_dens_curr_high, 1, mean)
str(pm.surv_JW_tmax_prev_dry_dens_curr_high)

# then calculate the credible intervals
CRI.surv_JW_tmax_prev_dry_dens_curr_high_95 <- apply(surv_JW_tmax_prev_dry_dens_curr_high, 1, function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_JW_tmax_prev_dry_dens_curr_high_95)
CRI.surv_JW_tmax_prev_dry_dens_curr_high_90 <- apply(surv_JW_tmax_prev_dry_dens_curr_high, 1, function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_JW_tmax_prev_dry_dens_curr_high_90)


surv_prob_JW_int2_tmax_prev_dry_dens_curr_high[[sex]] <- data.frame(state = "JW", 
                                                                    tmax_prev_dry_std = tmax_prev_dry_std,
                                                                    tmax_prev_dry = round(tmax_prev_dry, digits = 2),
                                                                    dens_curr_std = max_dens_std,
                                                                    dens_curr = round(max(dens)),
                                                                    dens_cat = "High Density",
                                                                    param = "survival",
                                                                    pred = pm.surv_JW_tmax_prev_dry_dens_curr_high,
                                                                    lower_95 = CRI.surv_JW_tmax_prev_dry_dens_curr_high_95[1,],
                                                                    upper_95 = CRI.surv_JW_tmax_prev_dry_dens_curr_high_95[2,],
                                                                    lower_90 = CRI.surv_JW_tmax_prev_dry_dens_curr_high_90[1,],
                                                                    upper_90 = CRI.surv_JW_tmax_prev_dry_dens_curr_high_90[2,],
                                                                    sex = sex)


surv_prob_JW_int2_tmax_prev_dry_dens_curr_high[[sex]]


### LOW DENSITY ###

# length.out rain_prev_wet X mcmc samples
surv_JW_tmax_prev_dry_dens_curr_low <- array(NA, dim = c(length(tmax.prev.stand), nrow(out)))

for (i in 1:length(tmax.prev.stand)){
  
  surv_JW_tmax_prev_dry_dens_curr_low[i,] <- plogis(mean.phi_JW 
                                                    + beta.JW.rain_current*0 # mean rain value = 0
                                                    + beta.JW.tmax_prev*tmax.prev.stand[i]
                                                    + beta.JW.dens_current*min_dens_std
                                                    + beta.JW.tmax_prev_dens_current*tmax.prev.stand[i]*min_dens_std)^7 # ^7 for seasonal pred of wet season
  
}

str(surv_JW_tmax_prev_dry_dens_curr_low)

# then we take the mean of the mcmc list
pm.surv_JW_tmax_prev_dry_dens_curr_low <- apply(surv_JW_tmax_prev_dry_dens_curr_low, 1, mean)
str(pm.surv_JW_tmax_prev_dry_dens_curr_low)

# then calculate the credible intervals
CRI.surv_JW_tmax_prev_dry_dens_curr_low_95 <- apply(surv_JW_tmax_prev_dry_dens_curr_low, 1, function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_JW_tmax_prev_dry_dens_curr_low_95)
CRI.surv_JW_tmax_prev_dry_dens_curr_low_90 <- apply(surv_JW_tmax_prev_dry_dens_curr_low, 1, function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_JW_tmax_prev_dry_dens_curr_low_90)


surv_prob_JW_int2_tmax_prev_dry_dens_curr_low[[sex]] <- data.frame(state = "JW", 
                                                                   tmax_prev_dry_std = tmax_prev_dry_std,
                                                                   tmax_prev_dry = round(tmax_prev_dry, digits = 2),
                                                                   dens_curr_std = min_dens_std,
                                                                   dens_curr = round(min(dens)),
                                                                   dens_cat = "Low Density",
                                                                   param = "survival",
                                                                   pred = pm.surv_JW_tmax_prev_dry_dens_curr_low,
                                                                   lower_95 = CRI.surv_JW_tmax_prev_dry_dens_curr_low_95[1,],
                                                                   upper_95 = CRI.surv_JW_tmax_prev_dry_dens_curr_low_95[2,],
                                                                   lower_90 = CRI.surv_JW_tmax_prev_dry_dens_curr_low_90[1,],
                                                                   upper_90 = CRI.surv_JW_tmax_prev_dry_dens_curr_low_90[2,],
                                                                   sex = sex)


surv_prob_JW_int2_tmax_prev_dry_dens_curr_low[[sex]]


surv_prob_JW_int2_tmax_prev_dry_dens_curr_cat[[sex]] <- rbind(surv_prob_JW_int2_tmax_prev_dry_dens_curr_low[[sex]], surv_prob_JW_int2_tmax_prev_dry_dens_curr_high[[sex]])

surv_JW_int2_plot <- ggplot(surv_prob_JW_int2_tmax_prev_dry_dens_curr_cat[[sex]])+
  geom_ribbon(aes(x= tmax_prev_dry, ymin = lower_90, ymax = upper_90), fill="blue",alpha=0.3) +
  geom_ribbon(aes(x= tmax_prev_dry, ymin = lower_95, ymax = upper_90), fill="lightblue",alpha=0.3) +
  geom_line(aes(tmax_prev_dry, pred), col="blue", lwd=3, linetype=1) +
  theme_classic()+
  facet_grid(~dens_cat) + 
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

surv_JW_int2_plot

# survival wet season: int + rainfall_current_Wet_season + tmax_prev_season + density_current + random year effect + tmax_prev_season:density_current

# AW ----------------------------------------------------------------------
mean.phi_AW <- out$mean.phi_AW
beta.AW.rain_current <- out$beta.AW.rain_current
beta.AW.tmax_prev <- out$beta.AW.tmax_prev
beta.AW.dens_current <- out$beta.AW.dens_current
beta.AW.tmax_prev_dens_current <- out$beta.AW.tmax_prev_dens_current

eps_phi_AW <- vector("list", length = length(unique(yearCov$yearCat)))
for (i in 1:length(unique(yearCov$yearCat))) {
  eps_phi_AW[[i]] <- out[,paste0("eps_phi_AW[", i, "]")]
}

###################################################################
## TMAX OF THE PREVIOUS DRY SEASON & POPDENS INTERACTION ----
###################################################################
# for the output data frames below
tmax_prev_dry_std = tmax.prev.stand
tmax_prev_dry = tmax.prev.stand*tmax_sd+tmax_mean

# length.out rain_prev_wet X length.oyt.dens_current X mcmc samples
surv_AW_tmax_prev_dry_dens_curr <- array(NA, dim = c(length(tmax.prev.stand), length(dens.stand), nrow(out)))

for (i in 1:length(tmax.prev.stand)){
  for (j in 1: length(dens.stand)){
    
    surv_AW_tmax_prev_dry_dens_curr[i,j,] <- plogis(mean.phi_AW 
                                                    + beta.AW.rain_current*0 # mean rain value = 0
                                                    + beta.AW.tmax_prev*tmax.prev.stand[i]
                                                    + beta.AW.dens_current*dens.stand[j]
                                                    + beta.AW.tmax_prev_dens_current*tmax.prev.stand[i]*dens.stand[j])^7 # ^7 for seasonal pred of wet season
  }
}

str(surv_AW_tmax_prev_dry_dens_curr)

# then we take the mean of the mcmc list
pm.surv_AW_tmax_prev_dry_dens_curr <- apply(surv_AW_tmax_prev_dry_dens_curr, c(1,2), mean)
str(pm.surv_AW_tmax_prev_dry_dens_curr)

# then calculate the credible intervals
CRI.surv_AW_tmax_prev_dry_dens_curr_95 <- apply(surv_AW_tmax_prev_dry_dens_curr, c(1,2), function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_AW_tmax_prev_dry_dens_curr_95)
CRI.surv_AW_tmax_prev_dry_dens_curr_90 <- apply(surv_AW_tmax_prev_dry_dens_curr, c(1,2), function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_AW_tmax_prev_dry_dens_curr_90)


surv_prob_AW_int2_tmax_prev_dry_dens_curr[[sex]] <- expand.grid(tmax_prev_dry_std = tmax_prev_dry_std, dens_curr_std = dens.stand)


surv_prob_AW_int2_tmax_prev_dry_dens_curr[[sex]]$tmax_prev_dry = rep(round(tmax_prev_dry, digits = 2), length = length(tmax.prev.stand))
surv_prob_AW_int2_tmax_prev_dry_dens_curr[[sex]]$dens_curr = rep(round(dens), each = length(tmax.prev.stand))

surv_prob_AW_int2_tmax_prev_dry_dens_curr[[sex]]$state = "AW"
surv_prob_AW_int2_tmax_prev_dry_dens_curr[[sex]]$param = "survival"
surv_prob_AW_int2_tmax_prev_dry_dens_curr[[sex]]$pred = c(pm.surv_AW_tmax_prev_dry_dens_curr)
surv_prob_AW_int2_tmax_prev_dry_dens_curr[[sex]]$lower_95 = c(CRI.surv_AW_tmax_prev_dry_dens_curr_95[1,,])
surv_prob_AW_int2_tmax_prev_dry_dens_curr[[sex]]$upper_95 = c(CRI.surv_AW_tmax_prev_dry_dens_curr_95[2,,])
surv_prob_AW_int2_tmax_prev_dry_dens_curr[[sex]]$lower_90 = c(CRI.surv_AW_tmax_prev_dry_dens_curr_90[1,,])
surv_prob_AW_int2_tmax_prev_dry_dens_curr[[sex]]$upper_90 = c(CRI.surv_AW_tmax_prev_dry_dens_curr_90[2,,])
surv_prob_AW_int2_tmax_prev_dry_dens_curr[[sex]]$sex = sex

surv_prob_AW_int2_tmax_prev_dry_dens_curr[[sex]]


AW_tmax_dens_int2_plot <- ggplot(surv_prob_AW_int2_tmax_prev_dry_dens_curr[[sex]])+
  geom_ribbon(aes(x= tmax_prev_dry, ymin = lower_90, ymax = upper_90, group = dens_curr , fill = dens_curr ), alpha=0.03) +
  geom_line(aes(tmax_prev_dry, pred, group = dens_curr, col = dens_curr ), lwd=3, linetype=1) +
  theme_classic()+
  scale_fill_continuous(type = "gradient", low = "#bcbddc", high = "#3f007d") +
  scale_color_continuous(type = "gradient", low = "#bcbddc", high = "#3f007d") +
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

AW_tmax_dens_int2_plot

##############################################################################
## TMAX OF THE PREVIOUS DRY SEASON & POPDENS INTERACTION CATEGORICAL ----
##############################################################################
# FOR HIGH DENSITY = 194 INDS, dens_std = 3.155485
# FOR LOW DENSITY = 40 INDS, dens_std = -1.322384

# for the output data frames below
tmax_prev_dry_std = tmax.prev.stand
tmax_prev_dry = tmax.prev.stand*tmax_sd+tmax_mean

### HIGH DENSITY ###

# length.out rain_prev_wet X mcmc samples
surv_AW_tmax_prev_dry_dens_curr_high <- array(NA, dim = c(length(tmax.prev.stand), nrow(out)))

for (i in 1:length(tmax.prev.stand)){
  
  surv_AW_tmax_prev_dry_dens_curr_high[i,] <- plogis(mean.phi_AW 
                                                     + beta.AW.rain_current*0 # mean rain value = 0
                                                     + beta.AW.tmax_prev*tmax.prev.stand[i]
                                                     + beta.AW.dens_current*max_dens_std
                                                     + beta.AW.tmax_prev_dens_current*tmax.prev.stand[i]*max_dens_std)^7 # ^7 for seasonal pred of wet season
  
}

str(surv_AW_tmax_prev_dry_dens_curr_high)

# then we take the mean of the mcmc list
pm.surv_AW_tmax_prev_dry_dens_curr_high <- apply(surv_AW_tmax_prev_dry_dens_curr_high, 1, mean)
str(pm.surv_AW_tmax_prev_dry_dens_curr_high)

# then calculate the credible intervals
CRI.surv_AW_tmax_prev_dry_dens_curr_high_95 <- apply(surv_AW_tmax_prev_dry_dens_curr_high, 1, function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_AW_tmax_prev_dry_dens_curr_high_95)
CRI.surv_AW_tmax_prev_dry_dens_curr_high_90 <- apply(surv_AW_tmax_prev_dry_dens_curr_high, 1, function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_AW_tmax_prev_dry_dens_curr_high_90)


surv_prob_AW_int2_tmax_prev_dry_dens_curr_high[[sex]] <- data.frame(state = "AW", 
                                                                    tmax_prev_dry_std = tmax_prev_dry_std,
                                                                    tmax_prev_dry = round(tmax_prev_dry, digits = 2),
                                                                    dens_curr_std = max_dens_std,
                                                                    dens_curr = round(max(dens)),
                                                                    dens_cat = "High Density",
                                                                    param = "survival",
                                                                    pred = pm.surv_AW_tmax_prev_dry_dens_curr_high,
                                                                    lower_95 = CRI.surv_AW_tmax_prev_dry_dens_curr_high_95[1,],
                                                                    upper_95 = CRI.surv_AW_tmax_prev_dry_dens_curr_high_95[2,],
                                                                    lower_90 = CRI.surv_AW_tmax_prev_dry_dens_curr_high_90[1,],
                                                                    upper_90 = CRI.surv_AW_tmax_prev_dry_dens_curr_high_90[2,],
                                                                    sex = sex)


surv_prob_AW_int2_tmax_prev_dry_dens_curr_high[[sex]]


### LOW DENSITY ###

# length.out rain_prev_wet X mcmc samples
surv_AW_tmax_prev_dry_dens_curr_low <- array(NA, dim = c(length(tmax.prev.stand), nrow(out)))

for (i in 1:length(tmax.prev.stand)){
  
  surv_AW_tmax_prev_dry_dens_curr_low[i,] <- plogis(mean.phi_AW 
                                                    + beta.AW.rain_current*0 # mean rain value = 0
                                                    + beta.AW.tmax_prev*tmax.prev.stand[i]
                                                    + beta.AW.dens_current*min_dens_std
                                                    + beta.AW.tmax_prev_dens_current*tmax.prev.stand[i]*min_dens_std)^7 # ^7 for seasonal pred of wet season
  
}

str(surv_AW_tmax_prev_dry_dens_curr_low)

# then we take the mean of the mcmc list
pm.surv_AW_tmax_prev_dry_dens_curr_low <- apply(surv_AW_tmax_prev_dry_dens_curr_low, 1, mean)
str(pm.surv_AW_tmax_prev_dry_dens_curr_low)

# then calculate the credible intervals
CRI.surv_AW_tmax_prev_dry_dens_curr_low_95 <- apply(surv_AW_tmax_prev_dry_dens_curr_low, 1, function(x) quantile(x, c(0.025, 0.975))) # 95%
str(CRI.surv_AW_tmax_prev_dry_dens_curr_low_95)
CRI.surv_AW_tmax_prev_dry_dens_curr_low_90 <- apply(surv_AW_tmax_prev_dry_dens_curr_low, 1, function(x) quantile(x, c(0.050, 0.950))) # 90%
str(CRI.surv_AW_tmax_prev_dry_dens_curr_low_90)


surv_prob_AW_int2_tmax_prev_dry_dens_curr_low[[sex]] <- data.frame(state = "AW", 
                                                                   tmax_prev_dry_std = tmax_prev_dry_std,
                                                                   tmax_prev_dry = round(tmax_prev_dry, digits = 2),
                                                                   dens_curr_std = min_dens_std,
                                                                   dens_curr = round(min(dens)),
                                                                   dens_cat = "Low Density",
                                                                   param = "survival",
                                                                   pred = pm.surv_AW_tmax_prev_dry_dens_curr_low,
                                                                   lower_95 = CRI.surv_AW_tmax_prev_dry_dens_curr_low_95[1,],
                                                                   upper_95 = CRI.surv_AW_tmax_prev_dry_dens_curr_low_95[2,],
                                                                   lower_90 = CRI.surv_AW_tmax_prev_dry_dens_curr_low_90[1,],
                                                                   upper_90 = CRI.surv_AW_tmax_prev_dry_dens_curr_low_90[2,],
                                                                   sex = sex)


surv_prob_AW_int2_tmax_prev_dry_dens_curr_low[[sex]]


surv_prob_AW_int2_tmax_prev_dry_dens_curr_cat[[sex]] <- rbind(surv_prob_AW_int2_tmax_prev_dry_dens_curr_low[[sex]], surv_prob_AW_int2_tmax_prev_dry_dens_curr_high[[sex]])

surv_AW_int2_plot <- ggplot(surv_prob_AW_int2_tmax_prev_dry_dens_curr_cat[[sex]])+
  geom_ribbon(aes(x= tmax_prev_dry, ymin = lower_90, ymax = upper_90), fill="blue",alpha=0.3) +
  geom_ribbon(aes(x= tmax_prev_dry, ymin = lower_95, ymax = upper_90), fill="lightblue",alpha=0.3) +
  geom_line(aes(tmax_prev_dry, pred), col="blue", lwd=3, linetype=1) +
  theme_classic()+
  facet_grid(~dens_cat) + 
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

surv_AW_int2_plot


###########################################################################
# 3. Save the outputs ----
###########################################################################

# RUN THIS AT THE END WHEN YOU ARE SURE YOU RAN BOTH FEMALE AND MALE RESULTS

# JD

str(surv_prob_JD_int2_tmax_curr_dry_dens_curr_cat)
surv_prob_JD_int2_tmax_curr_dry_dens_curr_cat_df <- do.call(rbind, surv_prob_JD_int2_tmax_curr_dry_dens_curr_cat)
rownames(surv_prob_JD_int2_tmax_curr_dry_dens_curr_cat_df) <- NULL
surv_prob_JD_int2_tmax_curr_dry_dens_curr_cat_df$sex_long <- ifelse(surv_prob_JD_int2_tmax_curr_dry_dens_curr_cat_df$sex == "F", "Female", "Male")
write.csv(surv_prob_JD_int2_tmax_curr_dry_dens_curr_cat_df,  file = paste0(outdir, "surv_prob_JD_int2_tmax_curr_dry_dens_curr_cat_df_seasonal_FM.csv"))

# AD

str(surv_prob_AD_int2_tmax_curr_dry_dens_curr_cat)
surv_prob_AD_int2_tmax_curr_dry_dens_curr_cat_df <- do.call(rbind, surv_prob_AD_int2_tmax_curr_dry_dens_curr_cat)
rownames(surv_prob_AD_int2_tmax_curr_dry_dens_curr_cat_df) <- NULL
surv_prob_AD_int2_tmax_curr_dry_dens_curr_cat_df$sex_long <- ifelse(surv_prob_AD_int2_tmax_curr_dry_dens_curr_cat_df$sex == "F", "Female", "Male")
write.csv(surv_prob_AD_int2_tmax_curr_dry_dens_curr_cat_df,  file = paste0(outdir, "surv_prob_AD_int2_tmax_curr_dry_dens_curr_cat_df_seasonal_FM.csv"))

# JW

str(surv_prob_JW_int2_tmax_prev_dry_dens_curr_cat)
surv_prob_JW_int2_tmax_prev_dry_dens_curr_cat_df <- do.call(rbind, surv_prob_JW_int2_tmax_prev_dry_dens_curr_cat)
rownames(surv_prob_JW_int2_tmax_prev_dry_dens_curr_cat_df) <- NULL
surv_prob_JW_int2_tmax_prev_dry_dens_curr_cat_df$sex_long <- ifelse(surv_prob_JW_int2_tmax_prev_dry_dens_curr_cat_df$sex == "F", "Female", "Male")
write.csv(surv_prob_JW_int2_tmax_prev_dry_dens_curr_cat_df,  file = paste0(outdir, "surv_prob_JW_int2_tmax_prev_dry_dens_curr_cat_df_seasonal_FM.csv"))

# AW
str(surv_prob_AW_int2_tmax_prev_dry_dens_curr_cat)
surv_prob_AW_int2_tmax_prev_dry_dens_curr_cat_df <- do.call(rbind, surv_prob_AW_int2_tmax_prev_dry_dens_curr_cat)
rownames(surv_prob_AW_int2_tmax_prev_dry_dens_curr_cat_df) <- NULL
surv_prob_AW_int2_tmax_prev_dry_dens_curr_cat_df$sex_long <- ifelse(surv_prob_AW_int2_tmax_prev_dry_dens_curr_cat_df$sex == "F", "Female", "Male")
write.csv(surv_prob_AW_int2_tmax_prev_dry_dens_curr_cat_df,  file = paste0(outdir, "surv_prob_AW_int2_tmax_prev_dry_dens_curr_cat_df_seasonal_FM.csv"))
