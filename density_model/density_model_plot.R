############################################################################
# Plot population abundance estimates from density model
# Creates mean_seasonal_density.csv used as density covariate in survival models
#
# Run this script after density_model_run.R has completed.
# Update the load() path below to point to your timestamped output file.
#
# Date: June 2026
# Authors: Dilsad Dagtekin and Dominik Behr
###########################################################################

## House keeping ----
##########################

rm(list = ls())
setwd("/path_to_files/density_model/")

## Load libraries ----
##########################

library(ggplot2)
library(lubridate)
library(MCMCvis)
library(coda)
library(boot)
library(viridis)
library(scales)
library(cowplot)


## Load model output ----
##########################

# Load the output file produced by density_model_run.R
# (update filename to match your timestamped output):
load("model_output/mlemur_posteriors_density.RData")
# Loaded object: run_js_model_null_popsize (coda MCMC list)

# Combine chains into one matrix:
df.post.dens <- do.call(rbind, lapply(run_js_model_null_popsize, as.matrix))


## Define plot theme ----
##########################

colBG     <- "white"
colTitles <- "black"
colLabels <- "#888888"
colGrid   <- "#eeeeee"
font      <- "Helvetica"
fontSize  <- 20

theme_general <- function(){
  theme_minimal() %+replace%
    theme(
      axis.text      = element_text(colour = colLabels, size = fontSize, family = font),
      axis.title.x   = element_text(colour = colTitles, size = fontSize, family = font,
                                    margin = margin(t = 15)),
      axis.title.y   = element_text(colour = colTitles, size = fontSize, family = font,
                                    margin = margin(r = 15), angle = 90),
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(colour = colGrid),
      axis.line   = element_blank(),
      axis.ticks  = element_blank(),
      plot.background = element_rect(fill = colBG, color = colBG, linewidth = 0),
      plot.title  = element_text(colour = colTitles, size = fontSize, family = font,
                                 hjust = 0, margin = margin(b = 15)),
      legend.position       = c(0.9, 0.9),
      legend.justification  = "right",
      legend.text  = element_text(colour = colLabels, size = fontSize, family = font),
      legend.title = element_text(colour = colTitles, size = fontSize, family = font),
      legend.key.width = unit(35, "pt")
    )
}


## Extract N[t] estimates ----
##########################

n_cols <- which(grepl(pattern = "N\\[", x = colnames(df.post.dens)))

df.inds.est <- data.frame(
  t              = seq_along(n_cols),
  n.inds.mean    = sapply(n_cols, function(X) mean(df.post.dens[, X])),
  n.inds.lower95 = sapply(n_cols, function(X) quantile(df.post.dens[, X], 0.025, na.rm = TRUE)),
  n.inds.upper95 = sapply(n_cols, function(X) quantile(df.post.dens[, X], 0.975, na.rm = TRUE)),
  n.inds.lower90 = sapply(n_cols, function(X) quantile(df.post.dens[, X], 0.050, na.rm = TRUE)),
  n.inds.upper90 = sapply(n_cols, function(X) quantile(df.post.dens[, X], 0.950, na.rm = TRUE)),
  n.inds.lower50 = sapply(n_cols, function(X) quantile(df.post.dens[, X], 0.250, na.rm = TRUE)),
  n.inds.upper50 = sapply(n_cols, function(X) quantile(df.post.dens[, X], 0.750, na.rm = TRUE))
)

# Assign dates: postwet occasions = May, prewet occasions = November
# First occasion is postwet 1994:
df.inds.est$date <- as.Date("1994-05-01")
for(i in 2:nrow(df.inds.est)){
  df.inds.est$date[i] <- floor_date(as.Date(6 * 31, df.inds.est$date[i - 1]), unit = "month")
}

head(df.inds.est)


## Save density estimates ----
##########################

write.csv(df.inds.est, file = "mean_seasonal_density.csv", row.names = FALSE)


## Plot abundance over time ----
##########################

colPlot <- "cornflowerblue"

p.N.ind <- ggplot() +
  geom_ribbon(data = df.inds.est,
              aes(x = date, ymin = n.inds.lower90, ymax = n.inds.upper90),
              fill = colPlot, alpha = 0.25) +
  geom_line(data  = df.inds.est,
            aes(x = date, y = n.inds.mean),
            color = colPlot, linetype = "solid", linewidth = 0.8) +
  geom_point(data = df.inds.est,
             aes(x = date, y = n.inds.mean),
             color = colPlot, fill = "white", size = 3, stroke = 1, shape = 21) +
  scale_y_continuous(name   = "Population abundance",
                     limits = c(0, 220), breaks = seq(0, 300, 50)) +
  scale_x_date(name   = NULL,
               limits = as.Date(c("1994-05-01", "2020-11-01"))) +
  theme_general()

p.N.ind

# Export as PDF:
ggsave(filename = "p.density.pdf",
       plot     = p.N.ind,
       width    = 10,
       height   = 5,
       family   = "Helvetica",
       device   = cairo_pdf,
       bg       = "transparent")
