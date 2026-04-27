#%%
#The goal in this script is to generate some diagnostic plots such as trace plots, estimate tables, R_eff, and Rhat
library(tidyverse)
library(ggplot2)
library(rstan)
library(Hmisc)
library(forestplot)
theme_set(theme_bw())
options(OutDec = "·")

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

#%%
io_set = "Main Results"
base_dir = getwd()

model_input_dir = file.path(base_dir, "Model Input", io_set)
model_output_dir = file.path(base_dir, "Model Output", io_set)
visuals_output_root_dir = file.path(base_dir, "Visuals Output", io_set)

#%%
model_name = "LogisticRegression"
suffix = ""
curr_sev_class_type = "1997type"

gen_diag_plots_tables = function(model_name, suffix, curr_sev_class_type, effects_type = "random", corr = TRUE){
    curr_output_filename = paste0("Results_", model_name, suffix, "_mean=0_sd=2_sd_mean=0.5_sdsd=2_", curr_sev_class_type, ".rds")
    curr_visuals_output_dir = file.path(visuals_output_root_dir, paste0(model_name, suffix), curr_sev_class_type)
    curr_model_fit = readRDS(file.path(model_output_dir, curr_output_filename))

    #%%
    #Generate Trace plots
    traceplot_output_dir = file.path(curr_visuals_output_dir, "DiagPlots")
    if (!dir.exists(traceplot_output_dir)) {
      dir.create(traceplot_output_dir, recursive = TRUE)
    }

    #Parameters to plot
    main_pars = c("beta_sero_prior", "beta_reg", "intercept")
	if(effects_type == "random"){
		main_pars = c(main_pars, "sigma")
	}
    all_pars = main_pars
	
    main_traceplot = traceplot(curr_model_fit, pars = main_pars) + 
		theme(legend.position = "top", legend.justification = "left")

	ggsave(main_traceplot, filename = file.path(traceplot_output_dir, paste0("Traceplot_MainPars_", curr_sev_class_type, ".png")), 
				width = 12, height = 8, units = "in", dpi = 300)
	
	if(effects_type == "random"){
		eps_pars = c("eps_vals")	
        all_pars = c(all_pars, eps_pars)
		eps_traceplot = traceplot(curr_model_fit, pars = eps_pars) + 
			theme(legend.position = "top", legend.justification = "left")
		
		ggsave(eps_traceplot, filename = file.path(traceplot_output_dir, paste0("Traceplot_Eps_", curr_sev_class_type, ".png")), 
			width = 16, height = 18, units = "in", dpi = 300)

		if(corr){
			rho_pars = c("Rho_1", "Rho_2")
            all_pars = c(all_pars, rho_pars)
			rho_traceplot = traceplot(curr_model_fit, pars = rho_pars)+
			theme(axis.text.x = element_text(angle = 90, size = 6), legend.position = "top", legend.justification = "left")
			ggsave(rho_traceplot, filename = file.path(traceplot_output_dir, paste0("Traceplot_Rho_", curr_sev_class_type, ".png")), 
						width = 16, height = 18, units = "in", dpi = 300)
		}

	}    
    #Generates tables with the parameter estimates, ESS, and Rhat values
    par_ests = summary(curr_model_fit, pars = all_pars)$summary %>% data.frame() %>%
                    select(mean, X2.5., X97.5., sd, se_mean, n_eff, Rhat)
    par_ests = par_ests %>% select(mean, X2.5., X97.5., sd, se_mean, n_eff, Rhat) %>% 
                    rename(Mean = mean, `2.5%` = X2.5., `97.5%` = X97.5., SD = sd, SE_Mean = se_mean, ESS = n_eff, Rhat = Rhat)
    write.csv(par_ests, file.path(curr_visuals_output_dir, paste0("Parameter_Estimates_", curr_sev_class_type, ".csv")), row.names = TRUE)

}

gen_diag_plots_tables("LogisticRegression", "", "1997type")
gen_diag_plots_tables("LogisticRegression", "", "2009type")
gen_diag_plots_tables("LogisticRegression", "", "hospitalisation")

gen_diag_plots_tables("LogisticRegression", "_no_unknown", "1997type")
gen_diag_plots_tables("LogisticRegression", "_no_unknown", "2009type")
gen_diag_plots_tables("LogisticRegression", "_no_unknown", "hospitalisation")

gen_diag_plots_tables("LogisticRegression", "_no_fried", "1997type")

gen_diag_plots_tables("LogisticRegression", "_no_narvaez", "1997type")

gen_diag_plots_tables("LogisticRegression", "_no_sabchareon", "1997type")


gen_diag_plots_tables("FE", "", "1997type", effects_type = "fixed")
gen_diag_plots_tables("FE", "", "2009type", effects_type = "fixed")
gen_diag_plots_tables("FE", "", "hospitalisation", effects_type = "fixed")

gen_diag_plots_tables("NoCorr", "", "1997type", corr = FALSE)
gen_diag_plots_tables("NoCorr", "", "2009type", corr = FALSE)
gen_diag_plots_tables("NoCorr", "", "hospitalisation", corr = FALSE)

