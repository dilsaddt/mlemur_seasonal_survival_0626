# mlemur_seasonal_survival_0626
Seasonal survival analysis of the gray mouse lemur population in Kirindy Forest, Madagascar.
Code for the multistate, Jolly-Seber capture-recapture analysis.
States are juveniles in dry season (JD), juveniles in wet season (JW), adults in dry season (AD), and adults in wet seasin (AW).
Males and females are analysed seperately to include sex effect on survival.
Environmental factors used in the model to check survival patterns under these factors are rainfall, maximum temperature.
Density-dependence also checked in the models, population density was obtained by a seperate analysis to account for the unceratinty in the population density over years of the study period.


## Files

custom_nimble_functions.R = custom functions that are needed to run nimble models

data/mlemur_example_data = data for running density models
data/dens_example; mlemur_example_data_females; mlemur_example_data_males = data for running additive and interaction models

*** Data here is representative just to make code running, not the actual data used in the analyses. ***

density_model/density_model_run.R = density model code
density_model/density_model_plot = plotting the density model output

model_run/mlemur_additive_github_0626.R = additive model code
model_run/mlemur_interaction1_github_0626.R = interaction model 1 code (rainfall:density)
model_run/mlemur_interaction2_github_0626.R = interaction model 2 code (maximum temperature:density)

model_diagnostics/mlemur_diagnostics_github_0626.R = model diagnostics for convergence and parameter distributions
here only diagnostics of additive model is shown as an example, but the process is the same for interaction models

prediction/mlemur_prediction_additive_github_0626.R = predictions and visualization from the additive model
prediction/mlemur_prediction_interaction1_github_0626.R = predictions and visualization from the interaction model 1 (rainfall:density)
prediction/mlemur_prediction_interaction1_github_0626.R = predictions and visualization from the interaction model 2 (maximum temperature:density)

## Software
R version  4.2.1 

NIMBLE version 0.13.2 (through R-package 'nimble')
