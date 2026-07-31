#99_ExtraPlots.R
#Code here is primarily for generating extra results for the adapted Peters' regression
#a linear regression where the dependent variable is epsilon/sigma and the predictor is 1/N.

#Imports ----
#%%
library(tidyverse)
library(ggplot2)
library(rstan)
library(Hmisc)
library(forestplot)
library(readxl)
library(metafor)
theme_set(theme_bw())

# Data Input ----
#%%
io_set = "Main Results"
base_dir = getwd()
setwd("..")
index_dir = getwd() #Where the review index is in
setwd(base_dir)


model_input_dir = file.path(base_dir, "Processed Data", io_set)
model_output_dir = file.path(base_dir, "Model Output", io_set)

index_df = read_excel(file.path(index_dir, "ReviewIndex_Final.xlsx"))


sev_class_types = c("1997type", "2009type", "hospitalisation")


visuals_output_dir = file.path(base_dir, "Visuals Output", io_set, "PublicationBias")
if(!dir.exists(visuals_output_dir)) {
  dir.create(visuals_output_dir, recursive = TRUE)
}
#%%
# curr_sev_class_type = sev_class_types[1]
curr_model_name = "LogisticRegression"
suffix = ""

for(curr_sev_class_type in sev_class_types) {
    curr_model_input = readRDS(file.path(model_input_dir, paste0("data_", curr_sev_class_type, suffix, ".rds")))
    curr_model_output = readRDS(file.path(model_output_dir, paste0("Results_", curr_model_name, "_mean=0_sd=2_sd_mean=0.5_sdsd=2_", curr_sev_class_type, suffix, ".rds")))

    #First, we generate a funnel plot with random effects on the x-axis (logit scale) and sample size on the y-axis with dashed line at 0
    theta_vals = summary(curr_model_output, pars = "theta")$summary %>% data.frame %>% rename(PointEst = mean, Lower = `X2.5.`, Upper = `X97.5.`) %>% select(PointEst, Lower, Upper)
    global_mean = theta_vals %>% pull(PointEst) %>% mean
    #global_mean = 0.5
    #Get the random effects from the eps_vals transformed parameter in Stan
    eps_vals = summary(curr_model_output, pars = "eps_vals")$summary %>% data.frame %>% rename(PointEst = mean, Lower = `X2.5.`, Upper = `X97.5.`) %>% select(PointEst, Lower, Upper)
    curr_outcome_df = curr_model_input$outcome_df %>% cbind(eps_vals)

    #%%
    # Peter's regression: here for each posterior sample, we fit a linear regression with random effects 
    # as dependent variable and 1/n as the independent variable

    #Get all the posterior samples for the random effects
    x_vals = 1/curr_outcome_df$N #The predictor is 1/N - the inverse sample size for each study
    eps_summary = summary(curr_model_output, pars = "eps_vals")$summary %>% data.frame %>% rename(PointEst = mean, Lower = `X2.5.`, Upper = `X97.5.`) %>% select(PointEst, Lower, Upper)
    sigma_summary = summary(curr_model_output, pars = "sigma")$summary %>% data.frame %>% rename(PointEst = mean, Lower = `X2.5.`, Upper = `X97.5.`) %>% select(PointEst, Lower, Upper)

    eps_samples = extract(curr_model_output, pars = "eps_vals")$eps_vals #We predict the size of the random effects using inverse sample size
    sigma_samples = extract(curr_model_output, pars = "sigma")$sigma #We predict the size of the random effects using inverse sample size

    #to do a weighted linear regression, we divide each epsilon by the appropriate sigma
    #We can do this by labelling each outcome with the index of the sigma we need (join outcome_df with scenario_df) 
    eps_matcher = curr_model_input$outcome_df
    curr_scenario_df = curr_model_input$scenario_df #ScenIndex will give the index of sigma needed for each epsilon
    eps_matcher = eps_matcher %>% left_join(curr_scenario_df %>% select(Scenario, ScenIndex), by = "Scenario") %>% mutate(EpsIndex = 1:nrow(.))

    eps_samples_weighted = eps_samples
    se_samples = matrix(NA, nrow = nrow(eps_samples), ncol = ncol(eps_samples))

    for(i in 1:nrow(eps_matcher)) {
        eps_ind = eps_matcher$EpsIndex[i]
        sigma_ind = eps_matcher$ScenIndex[i]
        eps_samples_weighted[, eps_ind] = eps_samples[, eps_ind] / sigma_samples[, sigma_ind]
        se_samples[, eps_ind] = sigma_samples[, sigma_ind]
    }

    #Regression per epsilon/sigma sample
    reg_results = apply(eps_samples_weighted, 1, function(eps_sample) {
        model_data = data.frame(eps_sample = eps_sample, x_vals = x_vals)
        model_fit = lm(eps_sample ~ x_vals, data = model_data)
        return(coef(model_fit))
    })

    res_mat = t(reg_results) %>% data.frame %>% rename(Intercept = `X.Intercept.`, Slope = x_vals)

    n_seq = seq(0, max(x_vals), length.out = 100)
    reg_lines = apply(res_mat, 1, function(coefs) {
        intercept = coefs["Intercept"]
        slope = coefs["Slope"]
        y_vals = intercept + slope * n_seq
        return(y_vals = y_vals)
    })
    shade_df = data.frame(n_seq = n_seq, 
                        low_shade = apply(reg_lines, 1, function(line) quantile(line, probs = 0.025)), 
                        up_shade = apply(reg_lines, 1, function(line) quantile(line, probs = 0.975)),
                        mean_shade = apply(reg_lines, 1, mean))

    regression_plot = curr_outcome_df %>% 
    ggplot(aes(x = 1/N)) +
        geom_pointrange(aes(y = PointEst, ymin = Lower, ymax = Upper)) + 
        geom_line(data = shade_df, aes(x = n_seq, y = mean_shade), color = "red") + 
        geom_ribbon(data = shade_df, aes(x = n_seq, ymin = low_shade, ymax = up_shade), fill = "red", alpha = 0.2) + 
        labs(title = paste("Reciprocal Sample Size Versus Random Effects", curr_sev_class_type), x = "1/Sample Size", y = "Random Effect (logit scale)") + 
        theme_bw()

    #Output summary stats from the regression
    reg_summary = res_mat %>% summarise(Intercept_Mean = mean(Intercept), Intercept_Lower = quantile(Intercept, probs = 0.025), Intercept_Upper = quantile(Intercept, probs = 0.975), 
                                        Slope_Mean = mean(Slope), Slope_Lower = quantile(Slope, probs = 0.025), Slope_Upper = quantile(Slope, probs = 0.975))


    # write.csv(reg_summary, file = file.path(visuals_output_dir, paste0("WeightedRegressionSummary_", curr_sev_class_type, suffix, ".csv")), row.names = FALSE)
    # ggsave(curr_funnel_plot, filename = file.path(visuals_output_dir, paste0("FunnelPlot_", curr_sev_class_type, suffix, ".png")), width = 8, height = 6)
    # ggsave(regression_plot, filename = file.path(visuals_output_dir, paste0("WeightedRegressionPlot_", curr_sev_class_type, suffix, ".png")), width = 8, height = 6)
}

