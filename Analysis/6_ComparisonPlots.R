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

visuals_output_dir = file.path(base_dir, "Visuals Output", io_set, "ComparisonPlots")
if(!dir.exists(visuals_output_dir)) {
    dir.create(visuals_output_dir, recursive = TRUE)
}

# I^2 Comparison plot ----
#%% I^2 Comparison 
#First, we create a plot comparing the I^2 values across models runs + with the different estimation methods

#Read in heterogeneity estimates from separate estimation
curr_model_name = "LogisticRegression"
data_suffix = ""

get_het_comp_df = function(curr_sev_class_type) {
    het_estimates = readRDS(file.path(het_dir, paste0("I2_Estimates_", curr_model_name, "_", curr_sev_class_type, ".rds")))

    #Read in model fit results to get the Stan-based estimates
    curr_model_fit = readRDS(file.path(model_output_dir, 
                    paste0("Results_", curr_model_name, "_mean=0_sd=2_sd_mean=0.5_sdsd=2_", curr_sev_class_type, ".rds")))
    curr_model_input = readRDS(file.path(model_input_dir, 
                    paste0("data_", curr_sev_class_type, data_suffix, ".rds")))

    curr_model_I2 = summary(curr_model_fit, pars = "I2")$summary %>% data.frame %>% mutate(ScenIndex = 1:nrow(.)) %>%
                        rename(PointEst = mean, Lower = X2.5., Upper = X97.5.)

    curr_model_I2 = curr_model_input$scenario_df %>% select(Scenario, ScenIndex) %>% left_join(curr_model_I2, by = "ScenIndex") %>% 
                        select(Scenario, PointEst, Lower, Upper) %>% mutate(Type = "Stan") %>% filter(!is.na(PointEst))

    curr_sep_I2 = het_estimates %>% filter(!is.na(I2_est)) %>% select(Scenario, I2_est) %>% rename(PointEst = I2_est) %>% mutate(Type = "Separate")

    curr_I2_ests = bind_rows(curr_model_I2, curr_sep_I2) %>% mutate(SevClass = curr_sev_class_type)
}
curr_model_I2_ests = lapply(sev_class_types, get_het_comp_df) %>% bind_rows()

    

#%%
I2_est_plot = curr_model_I2_ests %>% 
    ggplot(aes(x = Scenario, y = PointEst, ymin = Lower, ymax = Upper, color = Type)) +
    geom_pointrange(position = position_dodge(width = 0.5)) + 
    theme(element_text = element_text(size = 12), legend.position = "top", legend.justification = "left") + 
    coord_flip() +
    facet_wrap(~SevClass, scales = "free_y")

ggsave(I2_est_plot, filename = file.path(visuals_output_dir, "I2_Estimate_Comparison.png"), width = 12, height = 8)

# Pooled effects comparison plot ----
#Next we create a plot comparing the values of p, the pooled effect size
#%%

curr_model_name = model_names[[1]]
curr_sev_class_type = sev_class_types[[1]]
data_suffix = ""

get_pooled_effects = function(curr_model_name, curr_sev_class_type, data_suffix = "") {
    curr_model_fit = readRDS(file.path(model_output_dir, 
                    paste0("Results_", curr_model_name, "_mean=0_sd=2_sd_mean=0.5_sdsd=2_", curr_sev_class_type, ".rds")))


    scenario_df = readRDS(file.path(model_input_dir, paste0("data_", curr_sev_class_type, data_suffix, ".rds")))$scenario_df 
    #Get the pooled effect values in p 
    curr_p_summ = summary(curr_model_fit, pars = "p")$summary %>% data.frame %>% mutate(RegSeroPriorInd = 1:nrow(.)) %>%
                        rename(PointEst = mean, Lower = X2.5., Upper = X97.5.) #%>% left_join(scenario_df, by = "ScenIndex") %>% 
                        # select(Scenario, PointEst, Lower, Upper) %>% mutate(Model = curr_model_name)

    curr_p_summ = scenario_df %>% left_join(curr_p_summ %>% select(RegSeroPriorInd, PointEst, Lower, Upper), by = c("RegSeroPriorInd")) %>%
                        filter(!is.na(PointEst))

    to_ret = curr_p_summ %>% select(Scenario, PointEst, Lower, Upper) %>% mutate(Model = curr_model_name, SevClass = curr_sev_class_type)
    return(to_ret)         
}
# pooled_effects_df = expand.grid(Model = model_names, SevClass = sev_class_types) %>% 
#                         pmap_dfr(function(Model, SevClass) get_pooled_effects(Model, SevClass, data_suffix))

pooled_effects_params = expand.grid(Model = model_names, SevClass = sev_class_types) %>% mutate(data_suffix = ifelse(Model == "LogisticRegression_no_unknown", "_no_unknown", ""))

pooled_effects_df = pooled_effects_params %>% pmap_dfr(function(Model, SevClass, data_suffix) get_pooled_effects(Model, SevClass, data_suffix))

#%%
pooled_effects_plot = pooled_effects_df %>%
    ggplot(aes(x = Scenario, y = PointEst, ymin = Lower, ymax = Upper, color = Model)) +
    geom_pointrange(position = position_dodge(width = 0.7)) + 
    theme(element_text = element_text(size = 12), legend.position = "top", legend.justification = "left") + 
    coord_flip() +
    facet_wrap(~SevClass, scales = "free_y")

ggsave(pooled_effects_plot, filename = file.path(visuals_output_dir, "Pooled_Effects_Comparison.png"), width = 14, height = 8)
pooled_effects_plot


#%% SeroPrior Comparison Plot
