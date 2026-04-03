library(tidyverse)
library(rstan)
library(ggplot2)

options(scipen=999)
#options(mc.cores = parallel::detectCores())
options(mc.cores = 4)
rstan_options(auto_write = TRUE)
Sys.setenv(RSTUDIO = "1")


#This assumes Analysis is the active working directory
base_dir = getwd()
model_code_dir = file.path(base_dir, "Model Code")

io_set = "Main Results"
model_input_dir = file.path(base_dir, "Processed Data", io_set)
output_dir = file.path(base_dir, "Model Output", io_set)

#Create the output_dir if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#We re-run this code three times, once for each fo the severity classification systems listed below
sev_classes = c("1997type", "2009type", "hospitalisation")
model_names = c("LogisticRegressionModel", "FE_Model", "NoCorr_Model")


#Set these as global variables for now
coeff_prior_mean = 0
coeff_prior_sd = 2
sd_prior_mean = 0.5
sd_prior_sd = 2
n_tau_sim = 1000
fit_model = function(sev_class_type, model_name, data_suffix = "", iter = 10000, chains = 4, max_treedepth = 13, adapt_delta = 0.99) {
  #Read the data based on some inputted parameters
  data_filename = paste0("data_", sev_class_type, data_suffix, ".rds")
  input_data_path = file.path(model_input_dir, data_filename)
  input_data = readRDS(input_data_path)

  #Some other parameters added to data_list (currently hardcoded but could be made more flexible in the future)
  data_list = input_data$data_list
  data_list$coeff_prior_mean = coeff_prior_mean
  data_list$coeff_prior_sd = coeff_prior_sd

  data_list$sd_prior_mean = sd_prior_mean
  data_list$sd_prior_sd = sd_prior_sd

  data_list$n_tau_sim = n_tau_sim

  #Assume that the model code is in the model_code_dir and we just load based on the model_name input
  model = stan_model(file.path(model_code_dir, paste0(model_name, ".stan")))
  fit_obj = sampling(model, data = data_list, iter = iter, chains = chains, seed = 0,
                     control = list(max_treedepth = max_treedepth, adapt_delta = adapt_delta))

  return(fit_obj)
}

#test = fit_model("1997type", "LogisticRegressionModel", iter = 500, chains = 1)

orig_1997 = fit_model("1997type", "LogisticRegressionModel", iter = 10000, chains = 4)
orig_2009 = fit_model("2009type", "LogisticRegressionModel", iter = 10000, chains = 4)
orig_hospitalisation = fit_model("hospitalisation", "LogisticRegressionModel", iter = 10000, chains = 4)

fe_1997 = fit_model("1997type", "FE_Model", iter = 10000, chains = 4)
fe_2009 = fit_model("2009type", "FE_Model", iter = 10000, chains = 4)
fe_hospitalisation = fit_model("hospitalisation", "FE_Model", iter = 10000, chains = 4)

no_corr_1997 = fit_model("1997type", "NoCorr_Model", iter = 10000, chains = 4, max_treedepth = 14)
no_corr_2009 = fit_model("2009type", "NoCorr_Model", iter = 10000, chains = 4, max_treedepth = 14)
no_corr_hospitalisation = fit_model("hospitalisation", "NoCorr_Model", iter = 10000, chains = 4, max_treedepth = 14)

orig_1997_no_unknown = fit_model("1997type", "LogisticRegressionModel", iter = 10000, chains = 4, data_suffix = "_no_unknown")
orig_2009_no_unknown = fit_model("2009type", "LogisticRegressionModel", iter = 10000, chains = 4, data_suffix = "_no_unknown")
orig_hospitalisation_no_unknown = fit_model("hospitalisation", "LogisticRegressionModel", iter = 10000, chains = 4, data_suffix = "_no_unknown")


#Loop through groups of runs and output the fitted results to files 
orig_fits = c(orig_1997, orig_2009, orig_hospitalisation)
fe_fits = c(fe_1997, fe_2009, fe_hospitalisation)
no_corr_fits = c(no_corr_1997, no_corr_2009, no_corr_hospitalisation)
orig_no_unknown_fits = c(orig_1997_no_unknown, orig_2009_no_unknown, orig_hospitalisation_no_unknown)
for(i in 1:length(orig_fits)) {
  curr_fit = orig_fits[[i]]
  curr_fe_fit = fe_fits[[i]]
  curr_no_corr_fit = no_corr_fits[[i]]
  curr_no_unknown_fit = orig_no_unknown_fits[[i]]
  curr_sev_class_type = sev_classes[i]
  saveRDS(curr_fit, file.path(output_dir, paste0("Results_LogisticRegression_mean=", coeff_prior_mean,"_sd=", coeff_prior_sd, "_sd_mean=",  sd_prior_mean, "_sdsd=", sd_prior_sd, "_", curr_sev_class_type, ".rds")))
  saveRDS(curr_fe_fit, file.path(output_dir, paste0("Results_FE_mean=", coeff_prior_mean,"_sd=", coeff_prior_sd, "_sd_mean=",  sd_prior_mean, "_sdsd=", sd_prior_sd, "_", curr_sev_class_type, ".rds")))
  saveRDS(curr_no_corr_fit, file.path(output_dir, paste0("Results_NoCorr_mean=", coeff_prior_mean,"_sd=", coeff_prior_sd, "_sd_mean=",  sd_prior_mean, "_sdsd=", sd_prior_sd, "_", curr_sev_class_type, ".rds")))
  saveRDS(curr_no_unknown_fit, file.path(output_dir, paste0("Results_Orig_NoUnknown_mean=", coeff_prior_mean,"_sd=", coeff_prior_sd, "_sd_mean=",  sd_prior_mean, "_sdsd=", sd_prior_sd, "_", curr_sev_class_type, ".rds")))
}

#summ = summary(fit_model)$summary
#summ = summary(orig_1997)$summary
