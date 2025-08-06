library(tidyverse)
library(rstan)
library(ggplot2)
Sys.setenv(RSTUDIO = 1)

options(scipen=999)
options(mc.cores = parallel::detectCores())

base_dir = getwd()
model_code_dir = base_dir
model_input_dir = file.path(base_dir, "Model Inputs")
output_dir = file.path(base_dir, "Model Outputs")

data_filename = "cohort_input_no_min_cases.rds"

input_data_path = file.path(model_input_dir, data_filename)
input_data = readRDS(input_data_path)


data_list = input_data$data_list
data_list$coeff_prior_mean = 0
data_list$coeff_prior_sd = 2

data_list$sd_prior_mean = 0.5
data_list$sd_prior_sd = 2


model = stan_model(file.path(model_code_dir, "CohortLogisticRegression_FixedEffects.stan"))
fit_model = sampling(model, data = data_list, iter = 10000, chains = 4, seed = 0,
                     control = list(max_treedepth = 12, adapt_delta = 0.99))

saveRDS(fit_model, file.path(output_dir, paste0("FE_Results_LogisticRegression_mean=", data_list$coeff_prior_mean,"_sd=", data_list$coeff_prior_sd, "_sd_mean=",  data_list$sd_prior_mean, "_sdsd=", data_list$sd_prior_sd, "_min_cases=1.rds")))
summ = summary(fit_model)$summary
