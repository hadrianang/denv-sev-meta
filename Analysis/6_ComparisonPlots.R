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
data_suffixes = c("", "_no_unknown", "", "")
visuals_output_dir = file.path(base_dir, "Visuals Output", io_set, "ComparisonPlots")
if(!dir.exists(visuals_output_dir)) {
    dir.create(visuals_output_dir, recursive = TRUE)
}

# I^2 Comparison plot ----
#%% I^2 Comparison 
#First, we create a plot comparing the I^2 values across models runs + with the different estimation methods
#This is only for the main model

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

I2_est_plot = curr_model_I2_ests %>% 
    ggplot(aes(x = Scenario, y = PointEst, ymin = Lower, ymax = Upper, color = Type)) +
    geom_pointrange(position = position_dodge(width = 0.5)) + 
    theme(text = element_text(size = 12), legend.position = "top", legend.justification = "left") + 
    coord_flip() +
    facet_wrap(~SevClass, scales = "free_y")

ggsave(I2_est_plot, filename = file.path(visuals_output_dir, "I2_Estimate_Comparison.png"), width = 12, height = 8)

# Pooled effects comparison plot ----
#Next we create a plot comparing the values of p, the pooled effect size
#%%

get_pooled_effects = function(curr_model_name, curr_sev_class_type, data_suffix = "") {
    curr_model_fit = readRDS(file.path(model_output_dir, 
                    paste0("Results_", curr_model_name, data_suffix, "_mean=0_sd=2_sd_mean=0.5_sdsd=2_", curr_sev_class_type, ".rds")))


    scenario_df = readRDS(file.path(model_input_dir, paste0("data_", curr_sev_class_type, data_suffix, ".rds")))$scenario_df 
    #Get the pooled effect values in p 
    curr_p_summ = summary(curr_model_fit, pars = "p")$summary %>% data.frame %>% mutate(RegSeroPriorInd = 1:nrow(.)) %>%
                        rename(PointEst = mean, Lower = X2.5., Upper = X97.5.) #%>% left_join(scenario_df, by = "ScenIndex") %>% 
                        # select(Scenario, PointEst, Lower, Upper) %>% mutate(Model = curr_model_name)

    curr_p_summ = scenario_df %>% left_join(curr_p_summ %>% select(RegSeroPriorInd, PointEst, Lower, Upper), by = c("RegSeroPriorInd")) %>%
                        filter(!is.na(PointEst))

    to_ret = curr_p_summ %>% select(Scenario, PointEst, Lower, Upper) %>% mutate(Model = paste0(curr_model_name, data_suffix), SevClass = curr_sev_class_type)
    return(to_ret)         
}
# pooled_effects_df = expand.grid(Model = model_names, SevClass = sev_class_types) %>% 
#                         pmap_dfr(function(Model, SevClass) get_pooled_effects(Model, SevClass, data_suffix))

pooled_effects_params = expand.grid(Model = model_names, SevClass = sev_class_types, stringsAsFactors = FALSE) %>% mutate(data_suffix = ifelse(Model == "LogisticRegression_no_unknown", "_no_unknown", "")) %>% 
                            mutate(Model = ifelse(Model == "LogisticRegression_no_unknown", "LogisticRegression", Model))

pooled_effects_df = pooled_effects_params %>% pmap_dfr(function(Model, SevClass, data_suffix) get_pooled_effects(Model, SevClass, data_suffix))

#%%
pooled_effects_plot = pooled_effects_df %>%
    ggplot(aes(x = Scenario, y = PointEst, ymin = Lower, ymax = Upper, color = Model)) +
    geom_pointrange(position = position_dodge(width = 0.7)) + 
    theme(text = element_text(size = 12), legend.position = "top", legend.justification = "left") + 
    coord_flip() +
    facet_wrap(~SevClass, scales = "free_y")

ggsave(pooled_effects_plot, filename = file.path(visuals_output_dir, "Pooled_Effects_Comparison.png"), width = 14, height = 8)
pooled_effects_plot

#%% SeroPrior Comparison Plot

gen_or_comp_df = function(curr_model_name, curr_sev_class_type, data_suffix = "") {
    curr_model_input = readRDS(file.path(model_input_dir, paste0("data_", curr_sev_class_type, data_suffix, ".rds")))
    scenario_df = curr_model_input$scenario_df 
    char_mat_guide = curr_model_input$char_mat_guide
    curr_model_fit = readRDS(file.path(model_output_dir, 
                    paste0("Results_", curr_model_name, data_suffix, "_mean=0_sd=2_sd_mean=0.5_sdsd=2_", curr_sev_class_type, ".rds")))

    ref_region = curr_model_input$ref_region
    non_ref_region = ifelse(ref_region == "Asia", "Americas", "Asia")
    non_ref_seroprior = char_mat_guide %>% filter(CharMatIndex >0) %>% pull(SeroPriorExp) %>% as.character
    ref_serotype_exposure = curr_model_input$ref_serotype_exposure
    or_labels = c("Intercept", non_ref_seroprior, non_ref_region)


    process_or_beta = function(par_name){
        beta_or_est  = summary(curr_model_fit, pars = c(par_name))$summary %>% data.frame %>%
                    mutate(Label = or_labels) %>% select(Label, mean, X2.5., X97.5.) %>% 
                    rename(PointEst = mean, Lower = X2.5., Upper = X97.5.)
        return(beta_or_est)
    }
    log_beta_or_est = process_or_beta("log_or_beta") %>% select(-Label) %>% rename(LogPointEst = PointEst, LogLower = Lower, LogUpper = Upper)
    beta_or_est = process_or_beta("or_beta")

    beta_or_est = cbind(beta_or_est, log_beta_or_est)


    sero_prior_mat = summary(curr_model_fit, pars = "sero_prior_or_ratios")$summary %>% data.frame %>%
                    select(mean, X2.5., X97.5.) %>%
                    rename(PointEst = mean, Lower = X2.5., Upper = X97.5.) %>%
                    mutate(Name = rownames(.)) %>%
                    separate(Name, sep = "\\[|,|\\]", into = c(NA, "ROW", "COL", NA)) %>% 
                    mutate(ROW = as.numeric(ROW), COL = as.numeric(COL)) %>%
                    filter(COL > ROW) %>% 
                    left_join(char_mat_guide, by = c("ROW" = "CharMatIndex")) %>% #Join the indices based on the input char_mat_guide
                    rename(RowScenario = SeroPriorExp) %>% 
                    left_join(char_mat_guide, by = c("COL" = "CharMatIndex")) %>%
                    rename(ColScenario = SeroPriorExp) %>% 
                    select(-c(ROW, COL))

    log_sero_prior_mat  = summary(curr_model_fit, pars = "log_sero_prior_or_ratios")$summary %>% data.frame %>%
                    select(mean, X2.5., X97.5.) %>%
                    rename(LogPointEst = mean, LogLower = X2.5., LogUpper = X97.5.) %>% mutate(Name = rownames(.)) %>%
                    separate(Name, sep = "\\[|,|\\]", into = c(NA, "ROW", "COL", NA)) %>% 
                    mutate(ROW = as.numeric(ROW), COL = as.numeric(COL)) %>%
                    filter(COL > ROW) %>%
                    select(-c(ROW, COL))

    sero_prior_mat = cbind(sero_prior_mat, log_sero_prior_mat)

    versus_ref_portion = beta_or_est %>% filter(!Label %in% c("Intercept", non_ref_region)) %>%
                            mutate(ColScenario = ref_serotype_exposure) %>% rename(RowScenario = Label)

    #We add in the part of the matrix versus the reference class (the reference is the denominator of the OR)
    sero_prior_mat_vis = sero_prior_mat %>% rbind(versus_ref_portion) %>% 
                        mutate(Significance = ifelse(((Lower > 1) | (Upper < 1)),"Significant", "Not Significant")) %>% #Check for significant OR   
                        arrange(desc(Significance), desc(PointEst)) %>% 
                        mutate(DispVal = sprintf("%.2f [%.2f to %.2f]", PointEst, Lower, Upper),
                            LogDispVal = sprintf("%.2f [%.2f to %.2f]", LogPointEst, LogLower, LogUpper))

    to_ret = sero_prior_mat_vis %>% mutate(SevClass = curr_sev_class_type, Model = paste0(curr_model_name, data_suffix))
    return(to_ret)
}

or_params_df = expand.grid(Model = model_names, SevClass = sev_class_types, stringsAsFactors = FALSE) %>% 
                        mutate(data_suffix = ifelse(Model == "LogisticRegression_no_unknown", "_no_unknown", "")) %>%
                        mutate(Model = ifelse(Model == "LogisticRegression_no_unknown", "LogisticRegression", Model))
  
or_comparison_df = or_params_df %>% pmap_dfr(function(Model, SevClass, data_suffix) gen_or_comp_df(Model, SevClass, data_suffix))

#%%
gen_or_comp_plot = function(or_comparison_df, curr_sev_class_type, use_model_column = FALSE, main_name = "LogisticRegression") {
    or_to_vis_df = or_comparison_df
    #Optionally filter out the unknown scenarios
    or_to_vis_df = or_to_vis_df %>% filter(!grepl("Unknown", RowScenario) & !grepl("Unknown", ColScenario))
    #We can create a separate plot for each SevClass type
    curr_to_vis = or_to_vis_df %>% filter(SevClass == curr_sev_class_type) %>% mutate(Comparison = paste0(RowScenario, " vs. ", ColScenario))
    # Build the header row
    # Build table rows - one row per unique scenario combination
    table_data = curr_to_vis %>% 
        select(RowScenario, ColScenario, Model, Significance) %>%
        pivot_wider(names_from = Model, values_from = Significance) %>%
        arrange(desc(.data[[main_name]] == "Significant"), RowScenario, ColScenario)
  
    table_ordering = table_data %>% select(RowScenario, ColScenario) %>% mutate(OrderID = 1:nrow(.))

    # Get model names directly from the pivoted columns
    models = setdiff(colnames(table_data), c("RowScenario", "ColScenario"))
    # Build header from actual model names
    if (use_model_column) {
        model_labels = models
    }else{ #If we do not use the model column, we just use these hard coded values
        model_labels = c("Main", "Excl. Unknown", "Fixed Effects", "No Correlation")
    }
    header = c("Numerator", "Denominator", model_labels)

    # Build table text
    table_text <- rbind(
    header,
    as.matrix(cbind(
        as.character(table_data$RowScenario),
        as.character(table_data$ColScenario),
        sapply(models, function(m) as.character(table_data[[m]]))
    ))
    )
    # Build mean/lower/upper matrices - one column per Model, one row per scenario
    scenarios <- table_data %>% select(RowScenario, ColScenario)

    make_matrix = function(col) {
        curr_to_vis %>% 
            pivot_wider(
            id_cols = c(RowScenario, ColScenario),
            names_from = Model,
            values_from = {{ col }}
            ) %>%
            left_join(table_ordering, by = c("RowScenario", "ColScenario")) %>%
            arrange(OrderID) %>%
            select(all_of(models)) %>%
            as.matrix()
    }

    mean_mat  = make_matrix(PointEst)
    lower_mat = make_matrix(Lower)
    upper_mat = make_matrix(Upper)

    # Add NA header row (forestplot needs it to align with table header)
    mean_mat  = rbind(NA, mean_mat)
    lower_mat = rbind(NA, lower_mat)
    upper_mat = rbind(NA, upper_mat)

    # Plot
    ret_plot = forestplot(
        labeltext = table_text,
        mean = mean_mat,
        lower = lower_mat,
        upper = upper_mat,
        is.summary = c(TRUE, rep(FALSE, nrow(table_data))),
        col = fpColors(
        box = RColorBrewer::brewer.pal(length(models), "Set1"),
        line = RColorBrewer::brewer.pal(length(models), "Set1")
        ),
        clip = c(0, 5), 
        legend = models,
        xlab = "Odds Ratio",
        boxsize = 0.2,
        zero = 1, # For odds ratios, the null value is 1
    )|> fp_set_style(txt_gp = fpTxtGp(cex = 0.8, ticks = gpar(cex = 0.8), xlab = gpar(cex = 0.8))) |>
        fp_set_zebra_style("#EFEFEF")#  
    return(ret_plot)
}

save_plot = function(curr_plot, filepath, width, height, units, res, type = "png"){
    if(type == "png"){
	   png(filepath, width = width, height = height, units = units, res = res)
    }else if(type == "eps"){
	   #Defaults to inches as units
	   postscript(filepath, width = width, height = height)
    }else{
	   pdf(filepath, width = width, height = height)
    }
    print(curr_plot)
    dev.off()
}


or_plot_1997 = gen_or_comp_plot(or_comparison_df, "1997type")
or_plot_2009 = gen_or_comp_plot(or_comparison_df, "2009type")
or_plot_hosp = gen_or_comp_plot(or_comparison_df, "hospitalisation")

save_plot(or_plot_1997, file.path(visuals_output_dir, "OR_Comparison_1997.png"), width = 14, height = 8, units = "in", res = 300)
save_plot(or_plot_2009, file.path(visuals_output_dir, "OR_Comparison_2009.png"), width = 14, height = 8, units = "in", res = 300)
save_plot(or_plot_hosp, file.path(visuals_output_dir, "OR_Comparison_Hospitalised.png"), width = 14, height = 8, units = "in", res = 300)

#%%
#Plotting prediction intervals versus means
curr_model_name = model_names[[1]]
curr_sev_class_type = sev_class_types[[1]]
data_suffix = ""

gen_interval_comp_df = function(curr_model_name, curr_sev_class_type, data_suffix = "") {
    curr_model_fit = readRDS(file.path(model_output_dir, 
                        paste0("Results_", curr_model_name, "_mean=0_sd=2_sd_mean=0.5_sdsd=2_", curr_sev_class_type, ".rds")))


    scenario_df = readRDS(file.path(model_input_dir, paste0("data_", curr_sev_class_type, data_suffix, ".rds")))$scenario_df 
    #Get the pooled effect values in p 
    curr_p_summ = summary(curr_model_fit, pars = "p")$summary %>% data.frame %>% mutate(RegSeroPriorInd = 1:nrow(.)) %>%
                        rename(PointEst = mean, Lower = X2.5., Upper = X97.5.)
    curr_p_summ = scenario_df %>% left_join(curr_p_summ %>% select(RegSeroPriorInd, PointEst, Lower, Upper), by = c("RegSeroPriorInd")) %>%
                        filter(!is.na(PointEst))

    curr_pred_interval_summ = summary(curr_model_fit, pars = "pred_interval_p")$summary %>% data.frame %>% mutate(ScenIndex = 1:nrow(.)) %>%
                        rename(PointEst = mean, Lower = X2.5., Upper = X97.5.)
    curr_pred_interval_summ = scenario_df %>% left_join(curr_pred_interval_summ %>% select(ScenIndex, PointEst, Lower, Upper), by = c("ScenIndex")) %>%
                        filter(!is.na(PointEst))
    interval_merged = curr_p_summ %>% select(Scenario, PointEst, Lower, Upper) %>% mutate(Type = "Mean") %>% 
        rbind(curr_pred_interval_summ %>% select(Scenario, PointEst, Lower, Upper) %>% mutate(Type = "PredInterval")) %>% 
        mutate(SevClass = curr_sev_class_type, Model = curr_model_name)
    return(interval_merged)
}

rand_model_names = c("LogisticRegression", "LogisticRegression_no_unknown", "NoCorr")
interval_params_df = expand.grid(Model = rand_model_names, SevClass = sev_class_types) %>% 
                        mutate(data_suffix = ifelse(Model == "LogisticRegression_no_unknown", "_no_unknown", ""))
interval_comparison_df = interval_params_df %>% pmap_dfr(function(Model, SevClass, data_suffix) gen_interval_comp_df(Model, SevClass, data_suffix))

#%%
for(curr_model in rand_model_names) {
    curr_interval_comparison_df = interval_comparison_df %>% filter(Model == curr_model)
    curr_interval_plot = curr_interval_comparison_df %>%
        ggplot(aes(x = Scenario, y = PointEst, ymin = Lower, ymax = Upper, color = Type)) +
        geom_pointrange(position = position_dodge(width = 0.5)) + 
        theme(text = element_text(size = 12), legend.position = "top", legend.justification = "left") + 
        coord_flip() +
        facet_wrap(~SevClass, scales = "free_y") + ggtitle(paste0("Mean vs. Prediction Interval Comparison - ", curr_model))
    ggsave(curr_interval_plot, filename = file.path(visuals_output_dir, paste0("Interval_Comparison_", curr_model, ".png")), width = 14, height = 8)
}

#%%
#Comparing results when excluding influential studies
#First, comparing the pooled meta-analytic estimates 
curr_model_name = model_names[[1]]
curr_sev_class_type = sev_class_types[[1]]
data_suffix = "no_narvaez"
gen_influential_studies_p_df = function(curr_model_name, curr_sev_class_type, data_suffix = "") {
    curr_model_fit = readRDS(file.path(model_output_dir, 
                        paste0("Results_", curr_model_name, data_suffix, "_mean=0_sd=2_sd_mean=0.5_sdsd=2_", curr_sev_class_type, ".rds")))


    scenario_df = readRDS(file.path(model_input_dir, paste0("data_", curr_sev_class_type,  data_suffix, ".rds")))$scenario_df 
    #Get the pooled effect values in p 
    curr_p_summ = summary(curr_model_fit, pars = "p")$summary %>% data.frame %>% mutate(RegSeroPriorInd = 1:nrow(.)) %>%
                        rename(PointEst = mean, Lower = X2.5., Upper = X97.5.)
    curr_p_summ = scenario_df %>% left_join(curr_p_summ %>% select(RegSeroPriorInd, PointEst, Lower, Upper), by = c("RegSeroPriorInd")) %>%
                        filter(!is.na(PointEst))

    dataset_name = ifelse(data_suffix == "", "Full Dataset", data_suffix)
    curr_p_summ = curr_p_summ %>% select(Scenario, PointEst, Lower, Upper) %>% mutate(Model = curr_model_name, SevClass = curr_sev_class_type, DataSet = dataset_name)
    return(curr_p_summ)         
}

influential_studies_p_params = data.frame(data_suffix = c("", "_no_narvaez", "_no_sabchareon", "_no_fried", "_no_unknown")) %>% 
                                            mutate(SevClass = "1997type", Model = "LogisticRegression")
#influential_studies_p_params = influential_studies_p_params %>% rbind(data.frame(data_suffix = "", SevClass = "1997type", Model = "LogisticRegression_no_unknown"))
influential_studies_p_df = influential_studies_p_params %>% pmap_dfr(function(Model, SevClass, data_suffix) gen_influential_studies_p_df(Model, SevClass, data_suffix))

influential_studies_p_plot = influential_studies_p_df %>%
    ggplot(aes(x = Scenario, y = PointEst, ymin = Lower, ymax = Upper, color = DataSet)) +
    geom_pointrange(position = position_dodge(width = 0.7)) + 
    theme(text = element_text(size = 12), legend.position = "top", legend.justification = "left") + 
    coord_flip()
ggsave(influential_studies_p_plot, filename = file.path(visuals_output_dir, paste0("Influential_Studies_P_Comparison_", curr_model, ".png")), width = 10, height = 15)

#%%
influential_studies_or_df = influential_studies_p_params %>% pmap_dfr(function(Model, SevClass, data_suffix) gen_or_comp_df(Model, SevClass, data_suffix))
influential_studies_or_df  = influential_studies_or_df %>% mutate(Model = ifelse(Model == "LogisticRegression", "Main", 
                                        ifelse(Model == "LogisticRegression_no_narvaez", "Excl. Narvaez", 
                                            ifelse(Model == "LogisticRegression_no_sabchareon", "Excl. Sabchareon",
                                                ifelse(Model == "LogisticRegression_no_fried", "Excl. Fried",
                                                    ifelse(Model == "LogisticRegression_no_unknown", "Excl. Unknown", Model))))))

or_plot_influential = gen_or_comp_plot(influential_studies_or_df, "1997type", use_model_column = TRUE, main_name = "Main")


save_plot(or_plot_influential, file.path(visuals_output_dir, "OR_Comparison_Influential.png"), width = 16, height = 8, units = "in", res = 300)