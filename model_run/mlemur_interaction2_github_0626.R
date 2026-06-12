############################################################################
# Multi-state model with capture-recapture analysis
# with Jolly-Seber model framework
# Code adapted from Kery and Schaub, BPA Book, Chapter 10
# Part 10.3.2

# This script uses the "samper" data from capture-recapture surveys of the gray mouse lemur
# population in the Kirindy Forest, Madagascar, between 1994 and 2020.
# The aim of this script is to run the multistate JS models
# with 5 states
# Juvenile-Postwet = JD (1)
# Juvenile-Prewet  = JW (2)
# Adult-Postwet    = AD (3)
# Adult-Prewet     = AW (4)
# not seen         = NS (5)

# "Sample period" has two sampling periods "prewet" and "postwet" (short: samper).
# For samper capture history we decided that samplings are done before and after wet season.
# For that we decided that our "natural year" is from May to May.
# Dry season comes before wet season.
# Dry season is from May to Sep and wet season is from Oct to Apr.
# Samplings are done in:
# April-May-June (4-5-6)     = after wet season, entering dry season: postwet, category 1
# September-October-November (9-10-11) = before wet season, entering wet season: prewet, category 2
# postwet comes before prewet
# for normal year 1994, the order is postwet 1994, prewet 1994

# Transition from 1 to 2 (e.g., postwet 1994 to prewet 1994): survival of dry season of 1994 (phi^5)
# mid-point (May) to mid-point (October): 5 months

# Transition from 2 to 1 (e.g., prewet 1994 to postwet 1995): survival of wet season of 1994 (phi^7)
# mid-point (October) to mid-point (May): 7 months

# Interaction model II: maximum temperature x population density

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
  library(parallel)
  library(MCMCvis)
  library(ggplot2)
  library(lubridate)
}

load.libraries()

## Load data ----
##########################

sex <- "F" # F for females M for males

if(sex=="F"){
  # load FEMALE data:
  load("../data/mlemur_example_data_females.RData") # capture history and initial values
  
  samperCH_ms <- samperCH_ms_F
  samperCH_inits <- samperCH_inits_F
  
}else if(sex=="M"){
  # load MALE data:
  load("../data/mlemur_example_data_males.RData") # capture history and initial values
  
  samperCH_ms <- samperCH_ms_M
  samperCH_inits <- samperCH_inits_M
}

str(samperCH_ms) # capture history
str(samperCH_inits) # initial values for the latent state

# covariates:
str(rain_example)       # standardized rainfall (length = n.occasions + 1)
str(tmax_example)       # standardized maximum temperature (length = n.occasions + 1)
str(yearCov_example)    # year index for random effects

# density model output as density covariate
load("../data/dens_example.RData")
str(dens_example)       # population density from the density model
dens <- dens_example$n.inds.mean
# standardized vector of length n.occasions + 1
dens.stand <- (dens-mean(dens))/sd(dens) # standardize
round(mean(dens.stand),10)==0 & round(sd(dens.stand),10)==1

## Load custom NIMBLE functions ----
##########################

source("custom_nimble_functions.R")

## NIMBLE MODEL ----
##########################

js_model_interaction2 <- nimbleCode( {

  #--------------------------------------
  # Parameters:
  # phi: survival probability
  # gamma: removal entry probability
  # p: capture probability
  #--------------------------------------
  # States (S): 6 states
  # 1 not yet entered NE
  # 2 Juvenile Postwet JD
  # 3 Juvenile Prewet  JW
  # 4 Adult Postwet    AD
  # 5 Adult Prewet     AW
  # 6 dead
  #
  # Observations (O): 5 states
  # 1 seen as JD
  # 2 seen as JW
  # 3 seen as AD
  # 4 seen as AW
  # 5 not seen
  #
  #-------------------------------------------------
  # Covariates:
  # rain[t]   = previous wet season rainfall  -> dry season survival (JD, AD)
  # rain[t+1] = current wet season rainfall   -> wet season survival (JW, AW)
  # tmax[t]   = previous dry season tmax      -> wet season survival (JW, AW)
  # tmax[t+1] = current dry season tmax       -> dry season survival (JD, AD)
  # dens[t+1] = current population density    -> all states
  #-------------------------------------------------

  ### PRIORS AND CONSTRAINTS ###

  # entry probability
  gamma[1:(n.occasions-1)] ~ dunif_vector_inScalar(min = 0, max = 1, length = (n.occasions-1))

  # priors for intercepts of survival
  mean.phi_JD ~ dnorm(0, sd = 2)
  mean.phi_JW ~ dnorm(0, sd = 2)
  mean.phi_AD ~ dnorm(0, sd = 2)
  mean.phi_AW ~ dnorm(0, sd = 2)

  # priors for intercepts of detection
  mean.p_JD ~ dnorm(0, sd = 2)
  mean.p_JW ~ dnorm(0, sd = 2)
  mean.p_AD ~ dnorm(0, sd = 2)
  mean.p_AW ~ dnorm(0, sd = 2)

  # priors for beta coefficients of covariate effects on survival:
  # JD (dry season): rain of previous wet season, tmax of current dry season, density
  beta.JD.rain_prev    ~ dnorm(mean = 0, sd = 2)
  beta.JD.tmax_current ~ dnorm(mean = 0, sd = 2)
  beta.JD.dens_current ~ dnorm(mean = 0, sd = 2)

  # JW (wet season): rain of current wet season, tmax of previous dry season, density
  beta.JW.rain_current ~ dnorm(mean = 0, sd = 2)
  beta.JW.tmax_prev    ~ dnorm(mean = 0, sd = 2)
  beta.JW.dens_current ~ dnorm(mean = 0, sd = 2)

  # AD (dry season): rain of previous wet season, tmax of current dry season, density
  beta.AD.rain_prev    ~ dnorm(mean = 0, sd = 2)
  beta.AD.tmax_current ~ dnorm(mean = 0, sd = 2)
  beta.AD.dens_current ~ dnorm(mean = 0, sd = 2)

  # AW (wet season): rain of current wet season, tmax of previous dry season, density
  beta.AW.rain_current ~ dnorm(mean = 0, sd = 2)
  beta.AW.tmax_prev    ~ dnorm(mean = 0, sd = 2)
  beta.AW.dens_current ~ dnorm(mean = 0, sd = 2)

  # interaction terms: maximum temperature x population density
  beta.JD.tmax_current_dens_current ~ dnorm(mean = 0, sd = 2)
  beta.AD.tmax_current_dens_current ~ dnorm(mean = 0, sd = 2)
  beta.JW.tmax_prev_dens_current    ~ dnorm(mean = 0, sd = 2)
  beta.AW.tmax_prev_dens_current    ~ dnorm(mean = 0, sd = 2)

  # priors for random year effect standard deviations
  sigma_phi_JD ~ dunif(0, 5)
  sigma_phi_JW ~ dunif(0, 5)
  sigma_phi_AD ~ dunif(0, 5)
  sigma_phi_AW ~ dunif(0, 5)

  sigma_p_JD ~ dunif(0, 5)
  sigma_p_JW ~ dunif(0, 5)
  sigma_p_AD ~ dunif(0, 5)
  sigma_p_AW ~ dunif(0, 5)

  # random year effects (vectorized)
  eps_phi_JD[1:n.years] ~ dnorm_vector_inScalar(mean = 0, sd = sigma_phi_JD, length = n.years)
  eps_phi_JW[1:n.years] ~ dnorm_vector_inScalar(mean = 0, sd = sigma_phi_JW, length = n.years)
  eps_phi_AD[1:n.years] ~ dnorm_vector_inScalar(mean = 0, sd = sigma_phi_AD, length = n.years)
  eps_phi_AW[1:n.years] ~ dnorm_vector_inScalar(mean = 0, sd = sigma_phi_AW, length = n.years)

  eps_p_JD[1:n.years] ~ dnorm_vector_inScalar(mean = 0, sd = sigma_p_JD, length = n.years)
  eps_p_JW[1:n.years] ~ dnorm_vector_inScalar(mean = 0, sd = sigma_p_JW, length = n.years)
  eps_p_AD[1:n.years] ~ dnorm_vector_inScalar(mean = 0, sd = sigma_p_AD, length = n.years)
  eps_p_AW[1:n.years] ~ dnorm_vector_inScalar(mean = 0, sd = sigma_p_AW, length = n.years)

  # prior for starting proportion of JD (1) and AD (2) in year 1
  pi[1:2] ~ ddirch(alpha = alpha.pi[1:2])


  ### Define state-transition and observation matrices ###

  for (t in 1:(n.occasions-1)){

    # monthly survival, interaction model II (tmax x density):
    logit(phi_month_JD[t]) <- mean.phi_JD + beta.JD.rain_prev*rain[t] +
      beta.JD.tmax_current*tmax[t+1] +
      beta.JD.dens_current*dens[t+1] +
      beta.JD.tmax_current_dens_current*tmax[t+1]*dens[t+1] +
      eps_phi_JD[year[t]]

    logit(phi_month_JW[t]) <- mean.phi_JW + beta.JW.rain_current*rain[t+1] +
      beta.JW.tmax_prev*tmax[t] +
      beta.JW.dens_current*dens[t+1] +
      beta.JW.tmax_prev_dens_current*tmax[t]*dens[t+1] +
      eps_phi_JW[year[t]]

    logit(phi_month_AD[t]) <- mean.phi_AD + beta.AD.rain_prev*rain[t] +
      beta.AD.tmax_current*tmax[t+1] +
      beta.AD.dens_current*dens[t+1] +
      beta.AD.tmax_current_dens_current*tmax[t+1]*dens[t+1] +
      eps_phi_AD[year[t]]

    logit(phi_month_AW[t]) <- mean.phi_AW + beta.AW.rain_current*rain[t+1] +
      beta.AW.tmax_prev*tmax[t] +
      beta.AW.dens_current*dens[t+1] +
      beta.AW.tmax_prev_dens_current*tmax[t]*dens[t+1] +
      eps_phi_AW[year[t]]

    # seasonal survival (dry season = 5 months, wet season = 7 months):
    phi_JD[t] <- phi_month_JD[t]^5
    phi_AD[t] <- phi_month_AD[t]^5
    phi_JW[t] <- phi_month_JW[t]^7
    phi_AW[t] <- phi_month_AW[t]^7

    # capture null model (year random effect only):
    logit(p_JD[t]) <- mean.p_JD + eps_p_JD[year[t]]
    logit(p_JW[t]) <- mean.p_JW + eps_p_JW[year[t]]
    logit(p_AD[t]) <- mean.p_AD + eps_p_AD[year[t]]
    logit(p_AW[t]) <- mean.p_AW + eps_p_AW[year[t]]


    ## Define probabilities of state S(t+1) given S(t) ##

    # from NE state:
    ps[1,1:M,t,1] <- nimNumeric(length = M, value = 1-gamma[t])
    ps[1,1:M,t,2] <- nimNumeric(length = M, value = (equals(t,1)*pi[1] + step(t-2))*gamma[t]*dry_season[t])
    ps[1,1:M,t,3] <- nimNumeric(length = M, value = 0)
    ps[1,1:M,t,4] <- nimNumeric(length = M, value = equals(t,1)*pi[2]*gamma[t])
    ps[1,1:M,t,5] <- nimNumeric(length = M, value = 0)
    ps[1,1:M,t,6] <- nimNumeric(length = M, value = 0)

    # from JD state:
    ps[2,1:M,t,1] <- nimNumeric(length = M, value = 0)
    ps[2,1:M,t,2] <- nimNumeric(length = M, value = 0)
    ps[2,1:M,t,3] <- nimNumeric(length = M, value = phi_JD[t])
    ps[2,1:M,t,4] <- nimNumeric(length = M, value = 0)
    ps[2,1:M,t,5] <- nimNumeric(length = M, value = 0)
    ps[2,1:M,t,6] <- nimNumeric(length = M, value = 1-phi_JD[t])

    # from JW state:
    ps[3,1:M,t,1] <- nimNumeric(length = M, value = 0)
    ps[3,1:M,t,2] <- nimNumeric(length = M, value = 0)
    ps[3,1:M,t,3] <- nimNumeric(length = M, value = 0)
    ps[3,1:M,t,4] <- nimNumeric(length = M, value = phi_JW[t])
    ps[3,1:M,t,5] <- nimNumeric(length = M, value = 0)
    ps[3,1:M,t,6] <- nimNumeric(length = M, value = 1-phi_JW[t])

    # from AD state:
    ps[4,1:M,t,1] <- nimNumeric(length = M, value = 0)
    ps[4,1:M,t,2] <- nimNumeric(length = M, value = 0)
    ps[4,1:M,t,3] <- nimNumeric(length = M, value = 0)
    ps[4,1:M,t,4] <- nimNumeric(length = M, value = 0)
    ps[4,1:M,t,5] <- nimNumeric(length = M, value = phi_AD[t])
    ps[4,1:M,t,6] <- nimNumeric(length = M, value = 1-phi_AD[t])

    # from AW state:
    ps[5,1:M,t,1] <- nimNumeric(length = M, value = 0)
    ps[5,1:M,t,2] <- nimNumeric(length = M, value = 0)
    ps[5,1:M,t,3] <- nimNumeric(length = M, value = 0)
    ps[5,1:M,t,4] <- nimNumeric(length = M, value = phi_AW[t])
    ps[5,1:M,t,5] <- nimNumeric(length = M, value = 0)
    ps[5,1:M,t,6] <- nimNumeric(length = M, value = 1-phi_AW[t])

    # from dead state:
    ps[6,1:M,t,1] <- nimNumeric(length = M, value = 0)
    ps[6,1:M,t,2] <- nimNumeric(length = M, value = 0)
    ps[6,1:M,t,3] <- nimNumeric(length = M, value = 0)
    ps[6,1:M,t,4] <- nimNumeric(length = M, value = 0)
    ps[6,1:M,t,5] <- nimNumeric(length = M, value = 0)
    ps[6,1:M,t,6] <- nimNumeric(length = M, value = 1)


    ## Define probabilities of O(t) given S(t) ##

    po[1,1:M,t,1] <- nimNumeric(length = M, value = 0)
    po[1,1:M,t,2] <- nimNumeric(length = M, value = 0)
    po[1,1:M,t,3] <- nimNumeric(length = M, value = 0)
    po[1,1:M,t,4] <- nimNumeric(length = M, value = 0)
    po[1,1:M,t,5] <- nimNumeric(length = M, value = 1)

    po[2,1:M,t,1] <- nimNumeric(length = M, value = p_JD[t])
    po[2,1:M,t,2] <- nimNumeric(length = M, value = 0)
    po[2,1:M,t,3] <- nimNumeric(length = M, value = 0)
    po[2,1:M,t,4] <- nimNumeric(length = M, value = 0)
    po[2,1:M,t,5] <- nimNumeric(length = M, value = 1-p_JD[t])

    po[3,1:M,t,1] <- nimNumeric(length = M, value = 0)
    po[3,1:M,t,2] <- nimNumeric(length = M, value = p_JW[t])
    po[3,1:M,t,3] <- nimNumeric(length = M, value = 0)
    po[3,1:M,t,4] <- nimNumeric(length = M, value = 0)
    po[3,1:M,t,5] <- nimNumeric(length = M, value = 1-p_JW[t])

    po[4,1:M,t,1] <- nimNumeric(length = M, value = 0)
    po[4,1:M,t,2] <- nimNumeric(length = M, value = 0)
    po[4,1:M,t,3] <- nimNumeric(length = M, value = p_AD[t])
    po[4,1:M,t,4] <- nimNumeric(length = M, value = 0)
    po[4,1:M,t,5] <- nimNumeric(length = M, value = 1-p_AD[t])

    po[5,1:M,t,1] <- nimNumeric(length = M, value = 0)
    po[5,1:M,t,2] <- nimNumeric(length = M, value = 0)
    po[5,1:M,t,3] <- nimNumeric(length = M, value = 0)
    po[5,1:M,t,4] <- nimNumeric(length = M, value = p_AW[t])
    po[5,1:M,t,5] <- nimNumeric(length = M, value = 1-p_AW[t])

    po[6,1:M,t,1] <- nimNumeric(length = M, value = 0)
    po[6,1:M,t,2] <- nimNumeric(length = M, value = 0)
    po[6,1:M,t,3] <- nimNumeric(length = M, value = 0)
    po[6,1:M,t,4] <- nimNumeric(length = M, value = 0)
    po[6,1:M,t,5] <- nimNumeric(length = M, value = 1)

  } # t


  ### LIKELIHOOD ###

  z[1:M,1] <- 1

  for (t in 2:n.occasions){
    for (i in 1:M){
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1, 1:6])
      y[i,t] ~ dcat(po[z[i,t], i, t-1, 1:5])
    } # i
  } # t

})


# Prepare initial values for z (true state) ----
#############################

head(samperCH_inits)
str(samperCH_ms)
str(samperCH_inits)

my.z.init    <- cbind(rep(1, dim(samperCH_inits)[1]), samperCH_inits)
colnames(my.z.init)[1] <- "1993_2"
nz           <- dim(samperCH_ms)[1]
my.z.init.ms <- rbind(my.z.init, matrix(1, ncol = dim(my.z.init)[2], nrow = nz))

CH.du <- cbind(rep(5, dim(samperCH_ms)[1]), samperCH_ms)
colnames(CH.du)[1] <- "1993_2"
CH.ms <- CH.du

dry_season <- (1:(dim(CH.ms)[2]-1)) %% 2
alpha.pi   <- c(sum(CH.du[,2] == 1), sum(CH.du[,2] == 3)) / min(sum(CH.du[,2] == 1), sum(CH.du[,2] == 3))


# Bundle constants and data ----
#############################

str(nimble_constants <- list(
  n.occasions = dim(CH.ms)[2],
  M           = dim(CH.ms)[1],
  dry_season  = dry_season,
  n.years     = length(unique(yearCov_example$year)),
  year        = yearCov_example$yearCat,
  alpha.pi    = alpha.pi,
  rain        = rain_example,
  tmax        = tmax_example,
  dens        = dens.stand
))

str(nimble_data <- list(y = CH.ms))


# Initial values ----
#############################

seed <- as.integer(Sys.time())
set.seed(seed)

inits <- function(){list(
  mean.phi_JD = rnorm(1, 0, 1), mean.phi_JW = rnorm(1, 0, 1),
  mean.phi_AD = rnorm(1, 0, 1), mean.phi_AW = rnorm(1, 0, 1),
  mean.p_JD   = rnorm(1, 0, 1), mean.p_JW   = rnorm(1, 0, 1),
  mean.p_AD   = rnorm(1, 0, 1), mean.p_AW   = rnorm(1, 0, 1),

  beta.JD.rain_prev    = rnorm(1, 0, 1),
  beta.JD.tmax_current = rnorm(1, 0, 1),
  beta.JD.dens_current = rnorm(1, 0, 1),
  beta.JW.rain_current = rnorm(1, 0, 1),
  beta.JW.tmax_prev    = rnorm(1, 0, 1),
  beta.JW.dens_current = rnorm(1, 0, 1),
  beta.AD.rain_prev    = rnorm(1, 0, 1),
  beta.AD.tmax_current = rnorm(1, 0, 1),
  beta.AD.dens_current = rnorm(1, 0, 1),
  beta.AW.rain_current = rnorm(1, 0, 1),
  beta.AW.tmax_prev    = rnorm(1, 0, 1),
  beta.AW.dens_current = rnorm(1, 0, 1),

  beta.JD.tmax_current_dens_current = rnorm(1, 0, 1),
  beta.AD.tmax_current_dens_current = rnorm(1, 0, 1),
  beta.JW.tmax_prev_dens_current    = rnorm(1, 0, 1),
  beta.AW.tmax_prev_dens_current    = rnorm(1, 0, 1),

  sigma_phi_JD = runif(1, 0, 5), sigma_phi_JW = runif(1, 0, 5),
  sigma_phi_AD = runif(1, 0, 5), sigma_phi_AW = runif(1, 0, 5),
  sigma_p_JD   = runif(1, 0, 5), sigma_p_JW   = runif(1, 0, 5),
  sigma_p_AD   = runif(1, 0, 5), sigma_p_AW   = runif(1, 0, 5),

  eps_phi_JD = rnorm(nimble_constants$n.years, 0, 5),
  eps_phi_JW = rnorm(nimble_constants$n.years, 0, 5),
  eps_phi_AD = rnorm(nimble_constants$n.years, 0, 5),
  eps_phi_AW = rnorm(nimble_constants$n.years, 0, 5),
  eps_p_JD   = rnorm(nimble_constants$n.years, 0, 5),
  eps_p_JW   = rnorm(nimble_constants$n.years, 0, 5),
  eps_p_AD   = rnorm(nimble_constants$n.years, 0, 5),
  eps_p_AW   = rnorm(nimble_constants$n.years, 0, 5),

  gamma = runif(nimble_constants$n.occasions - 1, 0, 1),
  z     = my.z.init.ms
)}

params <- c(
  "mean.phi_JD", "mean.phi_JW", "mean.phi_AD", "mean.phi_AW",
  "mean.p_JD",   "mean.p_JW",   "mean.p_AD",   "mean.p_AW",
  "beta.JD.rain_prev",    "beta.JD.tmax_current", "beta.JD.dens_current",
  "beta.JW.rain_current", "beta.JW.tmax_prev",    "beta.JW.dens_current",
  "beta.AD.rain_prev",    "beta.AD.tmax_current", "beta.AD.dens_current",
  "beta.AW.rain_current", "beta.AW.tmax_prev",    "beta.AW.dens_current",
  "beta.JD.tmax_current_dens_current", "beta.AD.tmax_current_dens_current",
  "beta.JW.tmax_prev_dens_current",    "beta.AW.tmax_prev_dens_current",
  "gamma", "pi",
  "sigma_phi_JD", "sigma_phi_JW", "sigma_phi_AD", "sigma_phi_AW",
  "sigma_p_JD",   "sigma_p_JW",   "sigma_p_AD",   "sigma_p_AW",
  "phi_JD",       "phi_JW",       "phi_AD",        "phi_AW",
  "phi_month_JD", "phi_month_JW", "phi_month_AD",  "phi_month_AW",
  "eps_phi_JD",   "eps_phi_JW",   "eps_phi_AD",    "eps_phi_AW",
  "p_JD",         "p_JW",         "p_AD",          "p_AW",
  "eps_p_JD",     "eps_p_JW",     "eps_p_AD",      "eps_p_AW"
)


# Compile and run model ----
#############################

start1 <- Sys.time()

js_int2.model <- nimbleModel(code      = js_model_interaction2,
                             constants = nimble_constants,
                             data      = nimble_data,
                             inits     = inits(),
                             calculate = FALSE)

C.js_int2.model <- compileNimble(js_int2.model)
js_int2.config  <- configureMCMC(js_int2.model, monitors = params)
js_int2.MCMC    <- buildMCMC(conf = js_int2.config, useConjugacy = FALSE)
C.js_int2.MCMC  <- compileNimble(js_int2.MCMC, project = js_int2.model)

end1 <- Sys.time()
runTime.Compile <- end1 - start1
runTime.Compile

nchains  <- 4
niter    <- 50000
thin     <- 4
nburnin  <- 10000

seed <- as.integer(Sys.time())
set.seed(seed) # for reproducibility
seed.vec <- sapply(X = 1:nchains, FUN = function(X) seed + X - 1)

start2 <- Sys.time()

run_js_interaction2 <- runMCMC(mcmc             = C.js_int2.MCMC,
                               nchains          = nchains,
                               nburnin          = nburnin,
                               niter            = niter,
                               thin             = thin,
                               samplesAsCodaMCMC = TRUE,
                               setSeed          = seed.vec)

end2 <- Sys.time()
runTime.Sample <- end2 - start2
runTime.Sample

# Save output ----
#############################

output_name   <- "js_model_interaction2_tmax_dens"
timestamp.now <- gsub(x = gsub(x = as.character(now()), pattern = ":", replacement = ""), pattern = " ", replacement = "-")
nz <- dim(samperCH_ms)[1]
save(runTime.Compile, runTime.Sample, seed, run_js_interaction2,
     file = paste0("model_output/", timestamp.now, "_", output_name,
                   "_chains", nchains, "_iter", niter,
                   "_burnin", nburnin, "_thin", thin,
                   "_seed", paste(seed.vec, collapse = "+"),
                   "_naugment", nz, ".RData"))
