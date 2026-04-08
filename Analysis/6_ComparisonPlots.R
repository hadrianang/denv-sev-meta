#6_ComparisonPlots
#The goal with this code is to generate plots comparing the results of different model runs.

#Imports ----
#%%
library(tidyverse)
library(ggplot2)
library(rstan)
library(Hmisc)
library(forestplot)
library(readxl)
theme_set(theme_bw())

inv_logit = plogis
logit = qlogis


# Data Directories ----
#%%
base_dir = getwd()

io_set = "Main Results"

model_input_dir = file.path(base_dir, "Processed Data", io_set)
model_output_dir = file.path(base_dir, "Model Output", io_set)
het_dir = file.path(base_dir, "Heterogeneity Estimates", io_set)

sev_class_types =  c("1997type", "2009type", "hospitalisation")
model_names = c("LogisticRegression", "LogisticRegression_no_unknown", "FE", "NoCorr")


# I^2 Comparison plot ----
#%% I^2 Comparison 
#First, we create a plot comparing the I^2 values across models runs + with the different estimation methods

#Read in heterogeneity estimates from separate estimation
curr_model_name = "LogisticRegression"
curr_sev_class_type = "1997type"
het_estimates = readRDS(file.path(het_dir, paste0("I2_Estimates_", curr_model_name, "_", curr_sev_class_type, ".rds")))

#Read in model fit results to get the Stan-based estimates
curr_model_fit = readRDS(file.path(model_output_dir, 
                  paste0("Results_", curr_model_name, "_mean=0_sd=2_sd_mean=0.5_sdsd=2_", curr_sev_class_type, ".rds")))

curr_model_I2 = summary(curr_model_fit, pars = "I2")$summary