library(tidyverse)
library(rstan)
library(ggplot2)

options(scipen=999)
#options(mc.cores = parallel::detectCores())
options(mc.cores = 4)
#Sys.setenv(RSTUDIO = "1")

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
model_names = c(_LogisticRegressionModel", "FE_Model", "NoCorr_Model")



fit_model = function(sev_class_type, model_name, data_suffix = "", iter = 10000, chains = 4, max_treedepth = 12, adapt_delta = 0.99) {
  #Read the data based on some inputted parameters
  data_filename = paste0("data_", sev_class_type, data_suffix, ".rds")
  input_data_path = file.path(model_input_dir, data_filename)
  input_data = readRDS(input_data_path)

  #Some other parameters added to data_list (currently hardcoded but could be made more flexible in the future)
  data_list = input_data$data_list
  data_list$coeff_prior_mean = 0
  data_list$coeff_prior_sd = 2

  data_list$sd_prior_mean = 0.5
  data_list$sd_prior_sd = 2

  data_list$n_tau_sim = 1000

  #Assume that the model code is in the model_code_dir and we just load based on the model_name input
  model = stan_model(file.path(model_code_dir, paste0(model_name, ".stan")))
  fit_obj = sampling(model, data = data_list, iter = iter, chains = chains, seed = 0,
                     control = list(max_treedepth = max_treedepth, adapt_delta = adapt_delta))

  return(fit_obj)
}
#model = stan_model(file.path(model_code_dir, "LogisticRegressionModel.stan"))
#model = stan_model(file.path(model_code_dir, "FE_Model.stan"))

#model = stan_model(file.path(model_code_dir, "NoCorr_Model.stan"))
# fit_model = sampling(model, data = data_list, iter = 10000, chains = 4, seed = 0,
#                      control = list(max_treedepth = 12, adapt_delta = 0.99))

#Just for testing
fit_model = sampling(model, data = data_list, iter = 500, chains = 1, seed = 0,
                     control = list(max_treedepth = 12, adapt_delta = 0.99))

#saveRDS(fit_model, file.path(output_dir, paste0("Results_LogisticRegression_mean=", data_list$coeff_prior_mean,"_sd=", data_list$coeff_prior_sd, "_sd_mean=",  data_list$sd_prior_mean, "_sdsd=", data_list$sd_prior_sd, "_", sev_class_type, ".rds")))
summ = summary(fit_model)$summary

#summary(fit_model, pars = "p")$summary %>% View
