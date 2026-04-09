# 4. Severity Class Visuals
# The goal of this notebook is to generate visuals for a single severity class system (e.g. 1997-type). We assume that the model has already been fitted and heterogeneity values computed prior to running the code in this notebook. This notebook has to be run thrice, once for each of the severity classification systems in the study. See "5_Group Visuals.ipynb" for the visuals that combine results across all three analyses.

#Imports ----
# %%
library(tidyverse)
library(ggplot2)
library(rstan)
library(Hmisc)
library(forestplot)
theme_set(theme_bw())

inv_logit = plogis
logit = qlogis

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

save_plot_all_formats = function(curr_plot, dir_path, filename, width, height, units, res){
    png_path = file.path(dir_path, "png", paste0(filename, ".png"))
    eps_path = file.path(dir_path, "eps", paste0(filename, ".eps"))
    pdf_path = file.path(dir_path, "pdf", paste0(filename, ".pdf"))

    save_plot(curr_plot, png_path, width, height, units, res, type = "png")
    save_plot(curr_plot, eps_path, width, height, units, res, type = "eps")
    save_plot(curr_plot, pdf_path, width, height, units, res, type = "pdf")
}

# Paths and Inputs ----
# %%
base_dir = getwd()
io_set = "Main Results"
save_output = TRUE


model_input_dir = file.path(base_dir, "Processed Data", io_set)
model_output_dir = file.path(base_dir, "Model Output", io_set)
het_dir = file.path(base_dir, "Heterogeneity Estimates", io_set)




#Code should be run for each of the severity classifications
sev_class_types = c("1997type", "2009type", "hospitalisation")
model_names = c("LogisticRegression", "FE", "NoCorr")

sev_class_type = sev_class_types[[1]]
model_name = model_names[[1]]
res_name = "Orig"
effects_type = "random"
het_est = "stan"
suffix = ""
#%%
#sev_class_type = "2009type"
#sev_class_type = "hospitalisation"
generate_visuals = function(sev_class_type, model_name, suffix = "", effects_type = "random", het_est = "separate"){
    res_name = model_name #To simplify the naming of results, we set res_name as model_name
    visuals_output_dir = file.path(base_dir, "Visuals Output", io_set, res_name, sev_class_type)
    subgroup_forests_dir = file.path(visuals_output_dir, "ScenarioForests")

    input_data = file.path(model_input_dir, paste0("data_", sev_class_type, suffix, ".rds")) %>% readRDS
    model_results = file.path(model_output_dir, paste0("Results_", model_name, suffix, "_mean=0_sd=2_sd_mean=0.5_sdsd=2_", sev_class_type, ".rds")) %>% readRDS


    if(!dir.exists(visuals_output_dir)){
        dir.create(visuals_output_dir, recursive = TRUE)
        
        dir.create(file.path(visuals_output_dir, "png"))
        dir.create(file.path(visuals_output_dir, "eps"))
        dir.create(file.path(visuals_output_dir, "pdf"))

        dir.create(subgroup_forests_dir, recursive = TRUE)

        dir.create(file.path(subgroup_forests_dir, "png"))
        dir.create(file.path(subgroup_forests_dir, "eps"))
        dir.create(file.path(subgroup_forests_dir, "pdf"))
    }

    # Model Fit Visuals ----
    # %%
    # First, we generate a visual showing how the model fits to the data.
    #The theta values can be compared to the values in the input outcome_df
    outcome_df = input_data$outcome_df

    #First, we get our data and compute the 95% binomial confidence intervals
    fit_vis_df = outcome_df %>% select(CovidenceID, Scenario, N, Severe) %>%
                cbind(binconf(outcome_df$Severe, outcome_df$N, method = "exact", return.df = TRUE)) %>% 
                mutate(Type = "DATA") %>% select(-c(N, Severe))

    #We retrieve the model estimates, which are the values of theta from the model results
    theta_est = summary(model_results, pars = "theta")$summary %>% data.frame

    fit_vals = theta_est %>% select(mean, X2.5., X97.5.) %>% 
                rename(PointEst = mean, Lower = X2.5., Upper = X97.5.) %>% 
                cbind(fit_vis_df %>% select(CovidenceID, Scenario)) %>% 
                mutate(Type = "MODEL")

    fit_vis_df = rbind(fit_vis_df, fit_vals)

    #We rename the single study with no author for ease of visualisation (this was a WHO bulletin)
    fit_vis_df = fit_vis_df %>% mutate(CovidenceID = ifelse(CovidenceID == "#9", "#9 - WHO 1987", CovidenceID))
    fit_vis_df

    options(repr.plot.width = 25, repr.plot.height = 33)
    fit_vis_df %>% 
        ggplot(aes(x = reorder(CovidenceID, PointEst), y = PointEst, ymin = Lower, ymax = Upper, colour = Type)) +
            geom_pointrange(position = position_dodge(width = 0.6)) + 
            coord_flip() + 
            theme(text = element_text(size = 22), legend.position = "top", legend.justification = "left") +
            xlab("Study") + ylab("Effect Size") + ylim(0,1) + 
            facet_wrap(vars(Scenario), scales = "free", ncol = 4)

    if(save_output){
        ggsave(file.path(visuals_output_dir, "png", paste0("ModelFitPerStudy", ".png")), width = 25, height = 33, unit = "in", dpi = 800, limitsize = FALSE)
        ggsave(file.path(visuals_output_dir, "eps", paste0("ModelFitPerStudy", ".eps")), width = 25, height = 33, unit = "in", dpi = 800, limitsize = FALSE)
        ggsave(file.path(visuals_output_dir, "pdf", paste0("ModelFitPerStudy", ".pdf")), width = 25, height = 33, unit = "in", dpi = 800, limitsize = FALSE)
    }

    # Scenario Summary ----
    # Here we generate a visual showing the pooled estimate for each scenario.
    #%%
    ref_region = input_data$ref_region
    non_ref_region = ifelse(ref_region == "Asia", "Americas", "Asia")
    ref_serotype_exposure = input_data$ref_serotype_exposure
    ref_scenario = input_data$ref_scenario

    #Name of what cases are considered severe, used to label the forestplot column
    severe_name = "Severe"
    if(sev_class_type == "1997type"){
        severe_name = "DHF/DSS"
    }else if(sev_class_type == "2009type"){
        severe_name = "SD"
    }else{
        severe_name = "Hospitalised"
    }

    non_ref_seroprior = input_data$char_mat_guide %>% filter(CharMatIndex >0) %>% pull(SeroPriorExp) %>% as.character
    scenario_labels = c(ref_scenario, 
                        paste0(ref_region, "-", non_ref_seroprior),
                        paste0(non_ref_region, "-", c(ref_serotype_exposure, non_ref_seroprior))
                    )

    #Get the scenario estimates, which are given as p from the model
    p_estimates = summary(model_results, pars = "p")$summary %>% data.frame
    p_estimates = p_estimates %>% select(mean, X2.5., X97.5.) %>% 
        rename(PointEst = mean, Lower = X2.5., Upper = X97.5.) %>% mutate(Scenario = scenario_labels)

    scenario_df = input_data$scenario_df #%>% left_join(het_results %>% select(Scenario, I2_est), by = "Scenario")
    scenario_vis_df = scenario_df %>% left_join(p_estimates, by = "Scenario") %>% 
                        filter(!is.na(CharMatIndex)) %>% 
                        mutate(DispVal = sprintf("%.2f%% [%.2f to %.2f%%]", 100*PointEst, 100*Lower, 100*Upper)) #%>%
                        #mutate(DispHet = ifelse(is.na(I2_est), "-", sprintf("%.2f%%", I2_est * 100))) %>% 
                        
    labeltext_list = c("Scenario", "Pooled", "NumStudies", "N", "Severe", "DispVal")
    header_args = list(
                        Scenario = "Scenario", 
                        Pooled = "Pooling?",
                        NumStudies = "# Studies",
                        N = "N", 
                        Severe = severe_name, 
                        DispVal = "Est. Proportion [95% CrI]"
                    )

    if(effects_type == "random"){
        #If we are using a random effects model, retrieve the heterogeneity estimates 
        
        tau_estimates = summary(model_results, pars = "tau")$summary %>% data.frame %>% 
                        select(mean, X2.5., X97.5.) %>% rename(PointEst = mean, Lower = X2.5., Upper = X97.5.) %>%
                        mutate(ScenIndex = 1:nrow(.)) %>% mutate(DispTau = sprintf("%.2f [%.2f to %.2f]", PointEst, Lower, Upper))
    
        #Here we have the option of using the original heterogeneity estimates (which were computed separately from the model fitting) or the new heterogeneity estimates 
        #(which are computed from the model fit and may differ slightly from the original ones due to differences in how tau^2 is estimated in the model fit versus the original heterogeneity estimation process).
        if(het_est == "separate"){
            het_results = file.path(het_dir, paste0("I2_Estimates_", model_name, suffix,  "_", sev_class_type, ".rds")) %>% readRDS
            scenario_vis_df = scenario_vis_df %>% left_join(het_results %>% select(Scenario, I2_est), by = "Scenario") %>%
                        mutate(DispHet = ifelse(is.na(I2_est), "-", sprintf("%.2f%%", I2_est * 100))) 
        }else if(het_est == "stan"){
            het_results = summary(model_results, pars = "I2")$summary %>% data.frame %>% 
                        select(mean, X2.5., X97.5.) %>% mutate(ScenIndex = 1:nrow(.)) %>% mutate(DispHet = sprintf("%.2f%% [%.2f to %.2f%%]", mean * 100, X2.5. * 100, X97.5. * 100))
            scenario_vis_df = scenario_vis_df %>% left_join(het_results %>% select(ScenIndex, DispHet), by = "ScenIndex")
        }else{
            stop("Invalid option for het_est. Must be either 'separate' or 'stan'.")
        }
        
        scenario_vis_df = scenario_vis_df %>% left_join(tau_estimates %>% select(ScenIndex, DispTau), by = "ScenIndex") %>% select(-ScenIndex)
        labeltext_list = c(labeltext_list, "DispHet", "DispTau")
        header_args$DispHet = expression(I^{2})
        header_args$DispTau = expression(tau)
    }
    scenario_vis_df = scenario_vis_df %>% mutate(Pooled = ifelse(Inclusion == "Included", "Pooled", "Not Pooled")) 


    #%%
    label_df = scenario_vis_df %>% 
        arrange(Scenario) %>%
        select(all_of(labeltext_list))
    #Generate a forest plot specifically for this severity class
    options(repr.plot.width = 15, repr.plot.height = 13)
    scenario_plot = scenario_vis_df %>% arrange(Scenario) %>%
        forestplot(labeltext = label_df,
                    mean = PointEst, lower = Lower, upper = Upper,
                    xticks = seq(0, 1, by = 0.2), ci.vertices = TRUE, boxsize = 0.2,
                    align = c("l", "l", "l", "r", "r"), graphwidth = unit(1.8, "in")
                ) |>
        fp_add_lines(h_2 = gpar(lty = 2)) |>
        fp_set_style(box = "royalblue", line = gpar(col = "darkblue"), summary = "royalblue",
                    txt_gp = fpTxtGp(cex =1.2, ticks = gpar(cex = 1)))
    scenario_plot = do.call(fp_add_header, c(list(scenario_plot), header_args))
    scenario_plot = scenario_plot |>
        fp_decorate_graph(graph.pos = 7) |> fp_set_zebra_style("#EFEFEF")

    if(save_output){
        #save_plot(scenario_plot, file.path(visuals_output_dir, paste0("SubgroupForest_",  sev_class_type, ".png")), 15, 13, "in", 600)
        save_plot_all_formats(scenario_plot, visuals_output_dir, paste0("SubgroupForest_",  sev_class_type), 15, 13, "in", 600)
    }


    # Odds Ratios versus Reference ----
    # Here we generate a visual of the odds ratios for each effect versus the reference class.
    # %%
    or_labels = c("Intercept", non_ref_seroprior, non_ref_region)

    process_or_beta = function(par_name){
        beta_or_est  = summary(model_results, pars = c(par_name))$summary %>% data.frame %>%
                    mutate(Label = or_labels) %>% select(Label, mean, X2.5., X97.5.) %>% 
                    rename(PointEst = mean, Lower = X2.5., Upper = X97.5.)
        return(beta_or_est)
    }
    log_beta_or_est = process_or_beta("log_or_beta") %>% select(-Label) %>% rename(LogPointEst = PointEst, LogLower = Lower, LogUpper = Upper)
    beta_or_est = process_or_beta("or_beta")

    beta_or_est = cbind(beta_or_est, log_beta_or_est)

    #We don't need the intercept
    ref_seroexp_label = data.frame(Label = paste0("Versus ", ref_serotype_exposure), PointEst = NA, Lower = NA, Upper  = NA, LogPointEst = NA, LogLower = NA, LogUpper = NA, DispVal = NA, LogDispVal = NA)
    ref_region_label = data.frame(Label = paste0("Versus ", ref_region), PointEst = NA, Lower = NA, Upper  = NA, LogPointEst = NA, LogLower = NA, LogUpper = NA, DispVal = NA, LogDispVal = NA)
    beta_or_vis_df = beta_or_est %>% filter(Label != "Intercept") %>% 
                        mutate(DispVal = sprintf("%.2f [%.2f to %.2f]", PointEst, Lower, Upper)) %>% 
                        mutate(LogDispVal = sprintf("%.2f [%.2f to %.2f]", LogPointEst, LogLower, LogUpper))

    or_seroexp_part = beta_or_vis_df %>% filter(Label != non_ref_region) %>% mutate(Label = paste0("     ", Label))
    or_region_part = beta_or_vis_df %>% filter(Label == non_ref_region)%>% mutate(Label = paste0("     ", Label))

    beta_or_vis_df = rbind(ref_seroexp_label, or_seroexp_part, ref_region_label, or_region_part)

    #Add significance markers *
    beta_or_vis_df = beta_or_vis_df %>% mutate(DispVal = ifelse(((Lower > 1) | (Upper < 1)), paste0(DispVal, "*"), DispVal),
                                            LogDispVal = ifelse(((LogLower > 0) | (LogUpper < 0)), paste0(LogDispVal, "*"), LogDispVal))

    options(repr.plot.width = 10, repr.plot.height = 10)
    or_versus_ref_plot = beta_or_vis_df %>%forestplot(labeltext = c(Label, DispVal), zero = 0,
                    mean = LogPointEst, lower = LogLower, upper = LogUpper,
                    xticks = seq(-5, 5, by = 1), ci.vertices = TRUE, boxsize = 0.2,
                    align = c("l", "l"), graphwidth = unit(2.2, "in"), xlab = "Est. Log Odds Ratio [95% CrI]"
                ) |>
        fp_add_lines(h_2 = gpar(lty = 2)) |>
        fp_set_style(box = "royalblue", line = gpar(col = "darkblue"), summary = "royalblue", zero = "red",
                    txt_gp = fpTxtGp(cex =1.2, ticks = gpar(cex = 1), xlab = gpar(cex = 1.2))) |>
        fp_add_header(
                        Label = "Variable", 
                        DispVal = "Est. Odds Ratio [95% CrI]"
                        #LogDispVal = "Est. Log Odds Ratio [95% CrI]"
        ) |> fp_set_zebra_style("#EFEFEF")

    if(save_output){
        #save_plot(or_versus_ref_plot, file.path(visuals_output_dir, paste0("ORvsRef_",  sev_class_type, ".png")), 10, 10, "in", 600)
        save_plot_all_formats(or_versus_ref_plot, visuals_output_dir, paste0("ORvsRef_",  sev_class_type), 10, 10, "in", 600)
    }


    # Odds Ratio Pairwise Comparison (Serotype - Prior Exposure) ----
    # Here, we generate a visual showing pairs of serotype-prior exposure effects and the odds ratio for severity between them (including the reference class and whether the difference was significant or not).
    #%%
    char_mat_guide = input_data$char_mat_guide

    #We create a visual to show pairwise odds ratios (these would be in the sero_prior_or_ratios matrix)
    sero_prior_mat = summary(model_results, pars = "sero_prior_or_ratios")$summary %>% data.frame %>%
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

    log_sero_prior_mat  = summary(model_results, pars = "log_sero_prior_or_ratios")$summary %>% data.frame %>%
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
    sero_prior_mat_vis

    # %%
    options(repr.plot.width = 13, repr.plot.height = 18)
    axis_lim = 5
    sero_prior_comp_plot = sero_prior_mat_vis %>%
                    forestplot(labeltext = c(RowScenario, ColScenario, DispVal, Significance), zero = 0,
                    mean = LogPointEst, lower = LogLower, upper = LogUpper, clip = c(-axis_lim, axis_lim),
                    xticks = seq(-axis_lim, axis_lim, by = as.integer(axis_lim/5)), boxsize = 0.2,
                    align = c("l", "l", "r", "r"), graphwidth = unit(2, "in"), ci.vertices.height = 0.2, xlab = "Est. Log Odds Ratio [95% CrI]"
                ) |>
        fp_add_lines(h_2 = gpar(lty = 2)) |>
        fp_set_style(box = "royalblue", line = gpar(col = "darkblue"), zero = "red", summary = "royalblue",
                    txt_gp = fpTxtGp(cex =1.2, ticks = gpar(cex = 1), xlab = gpar(cex = 1.2))) |>
        fp_add_header(
                        RowScenario = "Numerator", 
                        ColScenario = "Denominator", 
                        DispVal = "Est. Odds Ratio [95% CrI]",
                        #LogDispVal = "Est. Log Odds Ratio [95% CrI]",
                        Significance = "Significance"
        ) |> fp_set_zebra_style("#EFEFEF") |> fp_decorate_graph(graph.pos = 4)
    if(save_output){
        #save_plot(sero_prior_comp_plot, file.path(visuals_output_dir, paste0("SeroPriorComp_",  sev_class_type, ".png")), 15, 18, "in", 600)
        save_plot_all_formats(sero_prior_comp_plot, visuals_output_dir, paste0("SeroPriorComp_",  sev_class_type), 15, 18, "in", 600)
    }
    # %%
    #Same plot as abive but without the unknowns
    options(repr.plot.width = 13, repr.plot.height = 12)
    axis_lim = 5
    sero_prior_comp_plot_known = sero_prior_mat_vis %>% filter(!(str_detect(RowScenario, "Unknown") | str_detect(ColScenario, "Unknown"))) %>%
                    forestplot(labeltext = c(RowScenario, ColScenario, DispVal, Significance), zero = 0, 
                    clip = c(-axis_lim, axis_lim), 
                    mean = LogPointEst, lower = LogLower, upper = LogUpper,
                    xticks = seq(-axis_lim, axis_lim, by = as.integer(axis_lim/5)), 
                    ci.vertices = FALSE, boxsize = 0.2, ci.vertices.height = 0.2, 
                    align = c("l", "l", "r", "r"), graphwidth = unit(2, "in"), xlab = "Est. Log Odds Ratio [95% CrI]"
                ) |>
        fp_add_lines(h_2 = gpar(lty = 2)) |>
        fp_set_style(box = "royalblue", line = gpar(col = "darkblue"), zero = "red", summary = "royalblue",
                    txt_gp = fpTxtGp(cex =1.2, ticks = gpar(cex = 1.1), xlab = gpar(cex = 1.3))) |>
        fp_add_header(
                        RowScenario = "Numerator", 
                        ColScenario = "Denominator", 
                        DispVal = "Est. Odds Ratio [95% CrI]",
                        #LogDispVal = "Est. Log Odds Ratio [95% CrI]",
                        Significance = "Significance"
        ) |> fp_set_zebra_style("#EFEFEF") |> fp_decorate_graph(graph.pos = 5)
    if(save_output){
        #save_plot(sero_prior_comp_plot_known, file.path(visuals_output_dir, paste0("SeroPriorComp_Known_",  sev_class_type, ".png")), 13, 12, "in", 600)
        save_plot_all_formats(sero_prior_comp_plot_known, visuals_output_dir, paste0("SeroPriorComp_Known_",  sev_class_type), 13, 12, "in", 600)
    }


    # Scenario Forests ----
    # In this section, we make individual forest plots for each of the scenarios where we pool data.
    # %%
    theta_df = summary(model_results, pars = "theta")$summary %>% data.frame %>% 
                    select(mean, X2.5., X97.5.) %>% rename(PointEst = mean, Lower = X2.5., Upper = X97.5.)

    outcome_ests = outcome_df %>% cbind(theta_df) %>% 
                    mutate(DispVal = sprintf("%.2f%% [%.2f to %.2f%%]", 100*PointEst, 100*Lower, 100*Upper))

    inc_scenarios = scenario_vis_df %>% filter(Inclusion != "Excluded")

    for(i in 1:nrow(inc_scenarios)){
        curr_row = inc_scenarios[i, ]
        curr_scen = curr_row$Scenario
        curr_scen_outcomes = outcome_ests %>% filter(Scenario == curr_scen)
        options(repr.plot.width = 13, repr.plot.height = nrow(curr_scen_outcomes))
        curr_scenario_plot = curr_scen_outcomes %>%
                        forestplot(labeltext = c(CovidenceID, Scenario, N, Severe, DispVal),
                        clip = c(-axis_lim, axis_lim), 
                        mean = PointEst, lower = Lower, upper = Upper,
                        xticks = seq(0, 1, by = 0.2), 
                        ci.vertices = FALSE, boxsize = 0.2, ci.vertices.height = 0.2, 
                        align = c("l", "l", "r", "r", "r"), graphwidth = unit(2, "in")
                    ) |>
            fp_add_lines(h_2 = gpar(lty = 2)) |>
            fp_set_style(box = "royalblue", line = gpar(col = "darkblue"), summary = "royalblue",
                        txt_gp = fpTxtGp(cex =1.2, ticks = gpar(cex = 1))) |>
            fp_add_header(
                            CovidenceID = "CovidenceID", 
                            Scenario = "Scenario", 
                            N = "N",
                            Severe = "Severe",
                            DispVal = "Est. Proportion [95% CrI]"
            ) |>  
            fp_append_row(
                    mean = curr_row %>% pull(PointEst),
                    lower = curr_row %>% pull(Lower),
                    upper = curr_row %>% pull(Upper),
                    CovidenceID = "Summary",
                    DispVal = curr_row %>% pull(DispVal),
                    is.summary = TRUE
            ) |>
            fp_set_zebra_style("#EFEFEF") |> fp_decorate_graph(graph.pos = 6)
        
        if(save_output){
            #save_plot(curr_scenario_plot, file.path(subgroup_forests_dir, paste0("Forest_", curr_scen, ".png")), 13, nrow(curr_scen_outcomes), "in", 600)
            save_plot_all_formats(curr_scenario_plot, subgroup_forests_dir, paste0("Forest_",  curr_scen), 13, nrow(curr_scen_outcomes), "in", 600)
        }
    }
}
#%%


generate_visuals("1997type", "LogisticRegression", het_est = "stan")
generate_visuals("2009type", "LogisticRegression", het_est = "stan")
generate_visuals("hospitalisation", "LogisticRegression", het_est = "stan")

#%%
generate_visuals("1997type", "NoCorr", het_est = "stan")
generate_visuals("2009type", "NoCorr", het_est = "stan")
generate_visuals("hospitalisation", "NoCorr", het_est = "stan")

#%%
generate_visuals("1997type", "FE", effects_type = "fixed")
generate_visuals("2009type", "FE", effects_type = "fixed")
generate_visuals("hospitalisation", "FE", effects_type = "fixed")