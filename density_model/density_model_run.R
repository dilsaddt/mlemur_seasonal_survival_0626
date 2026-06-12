############################################################################
# Null JS model to estimate population abundance (density covariate)
# with Jolly-Seber model framework
# Code adapted from Kery and Schaub, BPA Book, Chapter 10
#
# This script uses capture-recapture data of the gray mouse lemur
# population in Kirindy Forest, Madagascar, between 1994 and 2020.
#
# Model structure: Bernoulli JS model with data augmentation
# Binary capture history: 1 = captured (any state), 0 = not captured
# Season-specific survival (dry/postwet season = 1, wet/prewet season = 2)
#
# Output: MCMC posteriors of N[t] (population size per occasion)
# Run density_model_plot.R afterwards to create mean_seasonal_density.csv
#
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

## Load custom NIMBLE functions ----
##########################

source("../custom_nimble_functions.R")

## Load data ----
##########################

# The density model uses binary (0/1) capture history.
# Here we use an example data
# For your own data: replace with load("data/your_binary_ch_data.RData")

load("../data/mlemur_example_data.RData")

# capture data
str(samperCH_ms)
# the initial values
str(samperCH_inits)

## NIMBLE MODEL ----
##########################

js_model_null_popsize <- nimbleCode( {
  
  ### PRIORS AND CONSTRAINTS ###
  
  # priors for means
  mean.phi[1:2] ~ dnorm_vector_inScalar(mean = 0, sd = 2, length = 2)    # Prior for mean survival prob of the two seasons
  mean.p[1:2] ~ dnorm_vector_inScalar(mean = 0, sd = 2, length = 2)    # Prior for mean detection prob of the two seasons
  
  # priors for hyperparameters of random season effects
  sigma_phi[1:2] ~ dunif_vector_inScalar(min = 0, max = 5, length = 2)  # Prior for standard deviation of phi of the two seasons
  sigma_p[1:2] ~ dunif_vector_inScalar(min = 0, max = 5, length = 2)    # Prior for standard deviation of p of the two seasons
  
  # resulting prior for season-dependent survival and detection probabilities
  logit(p[1]) ~ dnorm(mean = mean.p[season[1]], sd = sigma_p[season[1]]) # detection probability
  for(t in 2:n.occasions){
    logit(phi[t-1]) ~ dnorm(mean = mean.phi[season[t]], sd = sigma_phi[season[t]]) # survival probability
    logit(p[t]) ~ dnorm(mean = mean.p[season[t]], sd = sigma_p[season[t]]) # detection probability
  }
  
  # prior for inclusion probability of individuals (i.e. data augmentation parameter for number of individuals)
  psi ~ dunif(min = 0, max = 1) 
  
  # Dirichlet prior for entry probabilities
  beta.b[1:n.occasions] ~ dgamma_vector_inScalar(shape = 1, rate = 1, length = n.occasions)
  
  # generate vector of conditional entry probabilities:
  nu[1:n.occasions] <- nimF.nu(beta.b = beta.b[1:n.occasions], n.occ = n.occasions) 
  
  
  ### LIKELIHOOD ###
  
  for(i in 1:M){ # we loop over all real and pseudo individuals
    
    ## state process ##
    
    w[i] ~ dbern(psi) # indicator variable, determining whether individual i is part of the population (w=1) or not (w=0)
    
    # first occasion:
    z[i,1] ~ dbern(nu[1]) # state of individual i at first occasion t=1; individuals is alive (z=1) or not (z=0)
    
    ## observation process ##
    y[i,1] ~ dbern(w[i]*z[i,1]*p[1]) # was individual i observed at occasion 1?
    
    # for all subsequent occasions:
    for(t in 2:n.occasions){
      
      # state process
      q[i,t-1] <- 1-z[i,t-1] # Availability for recruitment
      z[i,t] ~ dbern(phi[t-1]*z[i,t-1] + nu[t]*prod(q[i,1:(t-1)])) # true state at occasion t
      
      # observation process
      y[i,t] ~ dbern(prob = w[i]*z[i,t]*p[t])
      
    } # t-loop
    
    # estimate auxiliary parameters
    u[i,1:n.occasions] <- w[i]*z[i,1:n.occasions]
    entered[i,1:n.occasions] <- nimF.entered(n.occ = n.occasions, u = u[i,1:n.occasions])
    
    # for the estimation of superpopulation size
    Nind[i] <- sum(u[i,1:n.occasions]) # for how many occasions was individual i alive?
    Nalive[i] <- 1-equals(Nind[i], 0) # was individual i ever alive?
    
  } # M-loop
  
  for(t in 1:n.occasions){
    
    # derived population parameters:
    N[t] <- sum(u[1:M,t]) # actual population size
    E[t] <- sum(entered[1:M,t]) # number of entries
    
  } # t-loop
  
  # superpopulation size
  Nsuper <- sum(Nalive[1:M])
  
})

# if exists remove the artificial columns 1993_1 and 1993_2
dummy_cols <- c("1993_1", "1993_2")
if (all(dummy_cols %in% colnames(samperCH_ms))) {
  samperCH_ms    <- samperCH_ms[,    !colnames(samperCH_ms)    %in% dummy_cols]
  samperCH_inits <- samperCH_inits[, !colnames(samperCH_inits) %in% dummy_cols]
}

# prepare initial values for z (true state):
# for the simulated data we can use CH.sur directly
# for the real data I prepared scripts
# samperCH_inits_250325.R
head(samperCH_inits)
# n.occs <- ncol(samperCH_inits)
# # n.occs <- 10
# a <- which(rowMeans(samperCH_ms[, 1:n.occs])<5)
# random.inds <- sample(x = a, size = 100, replace = FALSE)
# samperCH_ms <- samperCH_ms[random.inds, 1:n.occs] # trial - shorter

str(samperCH_ms)

# samperCH_inits <- samperCH_inits[random.inds,1:n.occs] # trial - shorter

str(samperCH_inits)

# change capture history to 0's and 1's:
CH.ms <- samperCH_ms
CH.ms[] <- 0
CH.ms[which(samperCH_ms%in%c(1:4))] <- 1
table(CH.ms)

# change z.init to 0's and 1's:
my.z.init.ms <- samperCH_inits
my.z.init.ms[] <- 0
my.z.init.ms[which(samperCH_inits%in%c(2:5))] <- 1
table(my.z.init.ms)

# augment:
n.aug <- round(0.5*nrow(CH.ms))
# n.aug <- 1
CH.ms <- rbind(CH.ms, matrix(data = 0, nrow = n.aug, ncol = ncol(CH.ms)))

m.z.init.aug <- matrix(data = 0, nrow = n.aug, ncol = ncol(my.z.init.ms))
for(i in 1:nrow(m.z.init.aug)) m.z.init.aug[i,sample(x = 1:ncol(m.z.init.aug), size = 1)] <- 1
my.z.init.ms <- rbind(my.z.init.ms, my.z.init.ms)


dim(CH.ms)
dim(my.z.init.ms)

# season
season <- rep(c(1,2), dim(CH.ms)[2]-1)

# Bundle data
str(nimble_constants <- list(n.occasions = dim(CH.ms)[2]
                             , M = dim(CH.ms)[1]
                             , season = season))

str(nimble_data <- list(y = CH.ms))

inits <- function(){list(mean.phi = rnorm(2, 0.5, 0.5) 
                         , mean.p = rnorm(2, 0.5, 0.5)
                         , sigma_phi = runif(2, 0, 2)
                         , sigma_p = runif(2, 0, 2)
                         , beta.b = rgamma(shape = 1, rate = 1, n = nimble_constants$n.occasions)
                         , z = my.z.init.ms
                         , psi = runif(1, 0, 1)
)}

params <- c("mean.phi" ,"mean.p" 
            ,"sigma_phi" ,"sigma_p"
            , "psi"
            , "N", "E" , "Nsuper"
            , "nu", "beta.b"
            # , "z"
)

start1 <- Sys.time()

js_popsize.model <- nimbleModel(code = js_model_null_popsize
                                , constants = nimble_constants
                                , data = nimble_data
                                , inits = inits()
                                , calculate = FALSE # set to FALSE to reduce compilation time 
                                # --> disables the calculation of all deterministic 
                                #     nodes and log-likelihood 
)

C.js_popsize.model <- compileNimble(js_popsize.model)
js_popsize.config <- configureMCMC(js_popsize.model
                                   , monitors = params # which parameters to monitor
)

js_popsize.MCMC <- buildMCMC(conf = js_popsize.config
                             , useConjugacy = FALSE # to reduce compilation time --> disables the search for conjugate samplers 
)

C.js_popsize.MCMC <- compileNimble(js_popsize.MCMC, project = js_popsize.model)


end1 <- Sys.time()

runTime.Compile <- end1-start1
runTime.Compile


seed <- as.integer(Sys.time())
seed.vec <- sapply(X = 1:nchains, FUN = function(X) seed + X - 1) # for reproducibility

nchains <- 2
niter <- 10000 #30000 # 10000
thin <- 1
nburnin <- 0

start2 <- Sys.time()

run_js_model_null_popsize <- runMCMC(mcmc              = C.js_popsize.MCMC,
                                     nchains           = nchains,
                                     nburnin           = nburnin,
                                     niter             = niter,
                                     thin              = thin,
                                     samplesAsCodaMCMC = TRUE,
                                     setSeed           = seed.vec)

end2 <- Sys.time()

runTime.Sample <- end2-start2
runTime.Sample

# not putting any burn-in to see where it starts to converge
# not putting any thinning to get all samples basically


# save output of all chains:
output_name <- "js_model_null_popsize"
timestamp.now <- gsub(x = gsub(x = as.character(now()), pattern = ":", replacement = ""), pattern = " ", replacement = "-")
save(runTime.Compile, runTime.Sample, seed, run_js_model_null_popsize, file = paste0("results/",timestamp.now,"_",output_name,"_chains",nchains,"_iter",niter,"_burnin",nburnin,"_thin",thin,"_seed",seed,"_naugment",n.aug,".RData"))
timestamp.now


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# POSTERIORS AND CONVERGENCE ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# conditional entry probability
MCMCsummary(run_js_model_null_popsize, params = c("nu"), round=2)
MCMCtrace(run_js_model_null_popsize, params = c("nu"), pdf = FALSE, ind = TRUE, Rhat = TRUE, n.eff = TRUE, iter = dim(run_js_model_null_popsize[[1]])[1])

# Survival parameters
MCMCsummary(run_js_model_null_popsize, params = c("mean.phi","sigma_phi"), round=2)
MCMCtrace(run_js_model_null_popsize, params = c("mean.phi","sigma_phi"), pdf = FALSE, ind = TRUE, Rhat = TRUE, n.eff = TRUE, iter = dim(run_js_model_null_popsize[[1]])[1])

# Detection parameters
MCMCsummary(run_js_model_null_popsize, params = c("mean.p","sigma_p"), round=2)
MCMCtrace(run_js_model_null_popsize, params = c("mean.p","sigma_p"), pdf = FALSE, ind = TRUE, Rhat = TRUE, n.eff = TRUE, iter = dim(run_js_model_null_popsize[[1]])[1])

# Size of the super population: -> to check if we need more data augmentation
MCMCsummary(run_js_model_null_popsize, params = c("Nsuper","psi"), round=2)
MCMCtrace(run_js_model_null_popsize, params = c("Nsuper","psi"), pdf = FALSE, ind = TRUE, Rhat = TRUE, n.eff = TRUE, iter = dim(run_js_model_null_popsize[[1]])[1])


# Abundance estimates
MCMCsummary(run_js_model_null_popsize, params = c("N"), round=2)
MCMCsummary(run_js_model_null_popsize, params = c("E"), round=2)

