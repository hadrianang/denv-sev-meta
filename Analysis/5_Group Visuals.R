# 5. Group Visuals
# The goal of this notebook is to generate visuals involving the results of the three different severity class systems that are part of the analysis. This should be run after heterogeneity metrics are computed. 

# %%
library(tidyverse)
library(ggplot2)
library(rstan)
library(Hmisc)
library(forestplot)
library(readxl)
theme_set(theme_bw())
options(OutDec = "·")

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
  
	#If the file format directories do not exist, create them
	png_dir = file.path(dir_path, "png")
	eps_dir = file.path(dir_path, "eps")
	pdf_dir = file.path(dir_path, "pdf")
	if(!dir.exists(png_dir)){
		dir.create(png_dir, recursive = TRUE)
		dir.create(eps_dir, recursive = TRUE)
		dir.create(pdf_dir, recursive = TRUE)
	}
  
    png_path = file.path(png_dir, paste0(filename, ".png"))
    eps_path = file.path(eps_dir, paste0(filename, ".eps"))
    pdf_path = file.path(pdf_dir, paste0(filename, ".pdf"))

    save_plot(curr_plot, png_path, width, height, units, res, type = "png")
    save_plot(curr_plot, eps_path, width, height, units, res, type = "eps")
    save_plot(curr_plot, pdf_path, width, height, units, res, type = "pdf")
}

save_output = TRUE


# 1. Data Input ----
# %%
base_dir = getwd()
setwd("..")
index_dir = getwd() #Where the review index is in
setwd(base_dir)

# io_set = "Main Results"
io_set = "AddDeDupe"

model_input_dir = file.path(base_dir, "Processed Data", io_set)
model_output_dir = file.path(base_dir, "Model Output", io_set)
sev_class_types = c("1997type", "2009type", "hospitalisation")

# 2. Scenario Probability Visuals ----
# Functions to generate the pooled scenario visuals
# %%
# curr_input_data = model_inputs[[1]]
# curr_output_data = model_outputs[[1]]
# curr_sev_class_type = sev_class_types[[1]]
# #curr_het_results = het_results[[1]]
# curr_het_results = NULL
# effect_type = "random"

process_scenario_probs = function(curr_input_data, curr_output_data, curr_sev_class_type, curr_het_results, effect_type = "random"){
	#Probabilities are given in the estimates of p 
	ref_region = curr_input_data$ref_region
	non_ref_region = ifelse(ref_region == "Asia", "Americas", "Asia")
	ref_serotype_exposure = curr_input_data$ref_serotype_exposure
	ref_scenario = curr_input_data$ref_scenario
	non_ref_seroprior = curr_input_data$char_mat_guide %>% filter(CharMatIndex >0) %>% pull(SeroPriorExp) %>% as.character
	
	
	scenario_labels = c(ref_scenario, 
					paste0(ref_region, "-", non_ref_seroprior),
					paste0(non_ref_region, "-", c(ref_serotype_exposure, non_ref_seroprior))
					)

	#Get estimates of p
	p_estimates = summary(curr_output_data, pars = "p")$summary %>% data.frame %>% 
			select(mean, X2.5., X97.5.) %>% rename(PointEst = mean, Lower = X2.5., Upper = X97.5.) %>% 
			mutate(Scenario = scenario_labels)

	scenario_df = curr_input_data$scenario_df %>% select(ScenIndex, Scenario, Inclusion, NumStudies, Severe, N)

	if(effect_type == "random"){
		tau_estimates = summary(curr_output_data, pars = "tau")$summary %>% data.frame %>% 
					select(mean, X2.5., X97.5.) %>% rename(PointEst = mean, Lower = X2.5., Upper = X97.5.) %>%
					mutate(ScenIndex = 1:nrow(.)) %>% mutate(DispTau = sprintf("%.2f [%.2f to %.2f]", PointEst, Lower, Upper))
					

		#If curr_het_results is NULL then use the I^2 estimates from Stan rather than the separate ones loaded from files
		if(is.null(curr_het_results)){
			curr_het_results = summary(curr_output_data, pars = "I2")$summary %>% data.frame %>% 
						select(mean, X2.5., X97.5.) %>% mutate(ScenIndex = 1:nrow(.)) %>% mutate(DispHet = sprintf("%.2f%% [%.2f to %.2f%%]", mean * 100, X2.5. * 100, X97.5. * 100))
			
			#We compute the probability that the I^2 values are >= 75% (substantial heterogeneity per Cochrane guide)
			curr_het_draws = rstan::extract(curr_output_data, pars = "I2")$I2	
			prob_het_sub = colMeans(curr_het_draws >= 0.75)
			prob_het_sub_df = data.frame(ScenIndex = 1:length(prob_het_sub), ProbHetSub = prob_het_sub) %>%
											mutate(DispProbHetSub = sprintf("%.2f%%", ProbHetSub * 100))
			scenario_df = scenario_df %>% left_join(curr_het_results %>% select(ScenIndex, DispHet), by = "ScenIndex") %>% 
								left_join(prob_het_sub_df, by = "ScenIndex") %>% select(-ProbHetSub)
		}else{
			scenario_df = scenario_df %>% left_join(curr_het_results %>% select(Scenario, I2_est), by = "Scenario") %>%
						mutate(DispHet = ifelse(is.na(I2_est), "-", sprintf("%.2f%%", I2_est * 100))) %>% select(-I2_est)
		}
		#Join the tau estimates to the scenario_df then remove the ScenIndex column to simplify
		scenario_df = scenario_df %>% left_join(tau_estimates %>% select(ScenIndex, DispTau), by = "ScenIndex")
		scenario_df = scenario_df %>%
			mutate(
				DispHet = str_replace_all(DispHet, "\\.", "·"),
				DispTau = str_replace_all(DispTau, "\\.", "·"),
				DispProbHetSub = str_replace_all(DispProbHetSub, "\\.", "·")
			)
	}

	#Get scenarios and then join the p estimates to them, removing rows where we do not have an estimate
	scenario_df = scenario_df %>% 
			left_join(p_estimates, by = "Scenario") %>% mutate(SevClass = curr_sev_class_type) %>%
			filter(!is.na(PointEst)) %>% arrange(Scenario) %>% 
			mutate(DispVal = sprintf("%.2f%% [%.2f to %.2f%%]", 100*PointEst, 100*Lower, 100*Upper)) %>% select(-ScenIndex) %>%
			mutate(DispVal = str_replace_all(DispVal, "\\.", "·"), DispN = ifelse(N < 10000, as.character(N), str_replace(as.character(N), "(\\d)(?=(\\d{3})+$)", "\\1 ")))

	disp_sev_class = ""
	if(curr_sev_class_type == "1997type"){
		disp_sev_class = "1997-type"
	}else if(curr_sev_class_type == "2009type"){
		disp_sev_class = "2009-type"
	}else{
		disp_sev_class = "Hospitalisation"
	}
	#Create the header row
	header_row = data.frame(Scenario = disp_sev_class, Inclusion = NA,
					NumStudies = NA, Severe = NA, N = NA, DispN = NA, PointEst = NA, 
						Lower = NA, Upper = NA, DispVal = NA, #DispHet = NA, DispTau = NA,
					SevClass = "HEADER")
	if (effect_type == "random"){
		header_row = header_row %>% mutate(DispHet = NA, DispTau = NA, DispProbHetSub = NA)
	}
	scenario_df = rbind(header_row, scenario_df %>% mutate(Scenario = paste0("       ", Scenario))) 
	return(scenario_df)
}
# %%
gen_scenario_pooled_visual = function(model_inputs, model_outputs, sev_class_types, het_results, visuals_output_dir, effect_type = "random", save_output = save_output){
	#If het_results is NULL, we use the Stan estimated I^2 values
	if(is.null(het_results)){
       	scenario_probs = mapply(process_scenario_probs, model_inputs, model_outputs, sev_class_types, 
							MoreArgs = list(curr_het_results = NULL, effect_type = effect_type), SIMPLIFY = FALSE)
	}else{
		#If het_results is not NULL, use those loaded from a file separately
        scenario_probs = mapply(process_scenario_probs, model_inputs, model_outputs, sev_class_types, het_results, 
							MoreArgs = list(effect_type = effect_type), SIMPLIFY = FALSE)
	}

	merged_scenarios = do.call(rbind, scenario_probs)
	merged_forest_df = merged_scenarios
	merged_forest_df = merged_forest_df %>% mutate(Pooled = ifelse(Inclusion == "Included", "Pooled", 
									ifelse(Inclusion == "Excluded", "Not Pooled", NA))) %>% 
										mutate(N = ifelse(Inclusion == "Excluded", "-", N),
										DispN = ifelse(Inclusion == "Excluded", "-", DispN),
										Severe = ifelse(Inclusion == "Excluded", "-", Severe))
  
	labeltext_list = c("Scenario", "Pooled", "NumStudies", "DispN", "Severe", "DispVal")
    header_args = list(
                        Scenario = "Scenario", 
                        Pooled = "Pooling?",
                        NumStudies = "# Studies",
                        DispN = "N", 
                        Severe = "Severe", 
                        DispVal = "Est. Proportion [95% CrI]"
                    )
	if(effect_type == "random"){
		labeltext_list = c(labeltext_list, "DispHet", "DispTau", "DispProbHetSub")
		header_args$DispHet = expression(I^{2})
		header_args$DispTau = expression(tau)
		header_args$DispProbHetSub = "Prob. I² ≥ 75%"
	}

    label_df = merged_forest_df %>% 

        select(all_of(labeltext_list))
	
	options(repr.plot.width = 19, repr.plot.height = 21)
	merged_forest_plot = merged_forest_df %>% 
		forestplot(labeltext = label_df,
				mean = PointEst, lower = Lower, upper = Upper,
				xticks = seq(0, 1, by = 0.2), ci.vertices = TRUE, boxsize = 0.2,
				align = c("l", "l", "l", "r", "r", "l", "l"), graphwidth = unit(1.8, "in"), 
				) |>
		fp_add_lines(h_2 = gpar(lty = 2)) |>
		fp_set_style(box = "royalblue", line = gpar(col = "darkblue"), summary = "royalblue",
				txt_gp = fpTxtGp(cex =1.2, ticks = gpar(cex = 1)))
      
        merged_forest_plot = do.call(fp_add_header, c(list(merged_forest_plot), header_args))
   		merged_forest_plot = merged_forest_plot |>
        fp_decorate_graph(graph.pos = 7) |> fp_set_zebra_style("#EFEFEF")


	if(save_output){
		#save_plot(merged_forest_plot, file.path(visuals_output_dir, paste0("ScenarioForest.png")), 15, 21, "in", 600)
		save_plot_all_formats(merged_forest_plot, visuals_output_dir, "ScenarioForest", 19, 23, "in", 600)
	}
	#Return the merged forest df for use in the vaccine trial comparison visuals
	return(merged_forest_df)
}

# %% [markdown]
# Versus Reference OR Table ----
# Here we generate an odds ratio table versus the reference value. Note that this assumes the same reference class and region across all three severity class systems. These values are stored in the or_beta quantity from the model results.

# %%
#We retrieve the reference region and serotype-prior exposure from the first model input
#This assumes that this stays uniform throughout the analysis with the three severity class systems. 

process_or_vals = function(curr_input_data, curr_output_data, curr_sev_class_type, log = FALSE){
    ref_region = curr_input_data$ref_region
    non_ref_region = ifelse(ref_region == "Asia", "Americas", "Asia")
    ref_serotype_exposure = curr_input_data$ref_serotype_exposure
    ref_scenario = curr_input_data$ref_scenario
    non_ref_seroprior = curr_input_data$char_mat_guide %>% filter(CharMatIndex >0) %>% pull(SeroPriorExp) %>% as.character
    
    or_labels = c("Intercept", non_ref_seroprior, non_ref_region) #Labels for the odds ratios
    
    #Get te odds ratio estimates
    par_name = "or_beta"
    if(log){
	   par_name = "log_or_beta"
    }
    or_beta_est = summary(curr_output_data, pars = par_name)$summary %>% data.frame %>% 
				select(mean, X2.5., X97.5.) %>% 
				rename(PointEst = mean, Lower = X2.5., Upper = X97.5.) %>% 
				mutate(Label = or_labels, 
					  SevClass = curr_sev_class_type,
					  DispVal = sprintf("%.2f [%.2f to %.2f]", PointEst, Lower, Upper)) #Add in the label and the severity system
    
    return(or_beta_est)
}

gen_ref_or_table_visual = function(model_inputs, model_outputs, sev_class_types, visuals_output_dir, save_output = save_output){
	#This function assumes that all inputs have the same reference region and serotype-prior exposure, so we can just pull these from the first input data	
	curr_input_data = model_inputs[[1]]
  	ref_region = curr_input_data$ref_region
    non_ref_region = ifelse(ref_region == "Asia", "Americas", "Asia")
    ref_serotype_exposure = curr_input_data$ref_serotype_exposure
    ref_scenario = curr_input_data$ref_scenario
  
  	sev_class_ors = mapply(process_or_vals, model_inputs, model_outputs, sev_class_types, SIMPLIFY = FALSE)
	log_sev_class_ors = mapply(process_or_vals, model_inputs, model_outputs, sev_class_types, rep(TRUE, 3), SIMPLIFY = FALSE)

	options(repr.plot.width = 13, repr.plot.height = 8)
	or_vis_df = do.call(rbind, sev_class_ors) %>% 
			mutate(DispVal = ifelse(((Lower > 1) | (Upper < 1) & Label != "Intercept"), paste0("***", DispVal), DispVal)) %>% #Check significance
			select(Label, SevClass, DispVal) %>% 
			pivot_wider(id_cols = Label, names_from = SevClass, values_from = DispVal)

	disp_order = expand.grid(c(paste0("DENV", 1:4)), c("Unknown", "Primary", "Secondary")) %>% 
		mutate(Order = paste0(Var2, "-", Var1)) %>% filter(Order != ref_serotype_exposure) %>% pull(Order)
	disp_order = c("Intercept", disp_order, non_ref_region)
	or_vis_df = or_vis_df %>% mutate(Label = factor(Label, levels = disp_order)) %>% arrange(Label)

	or_vis_df %>% mutate(mean = 0, lower = 0, upper = 0) %>%
		forestplot(labeltext = c(Label, `1997type`, `2009type`, hospitalisation), 
				xticks = seq(0, 1, by = 0.2), ci.vertices = TRUE, boxsize = 0.0,
				align = c("l", "r", "r", "r"),  graphwidth = unit(0, "cm")
				) |>
		fp_add_lines(h_2 = gpar(lty = 2)) |>
		fp_set_style(box = "royalblue", line = gpar(col = "darkblue"), summary = "royalblue",
				txt_gp = fpTxtGp(cex =1.2, ticks = gpar(cex = 0))) |>
		fp_add_header(
					Label = "Name",
					`1997type` = "1997-type",
					`2009type` = "2009-type",
					hospitalisation = "Hospitalisation") |> 
		#fp_decorate_graph(graph.pos = 0) |>
		fp_set_zebra_style("#EFEFEF")

	# %%
	or_vis_alt_df = do.call(rbind, sev_class_ors) %>% 
				pivot_wider(id_cols = c(Label), 
						names_from = SevClass, 
						values_from = c(PointEst, Lower, Upper, DispVal)) %>%
				mutate(Label = paste0("    ", Label))


	log_or_vis_alt_df = do.call(rbind, log_sev_class_ors) %>% 
				pivot_wider(id_cols = c(Label), 
						names_from = SevClass, 
						values_from = c(PointEst, Lower, Upper, DispVal))%>%
				mutate(Label = paste0("    ", Label))


	ref_seroprior_label = as.data.frame(as.list(rep(NA, ncol(or_vis_alt_df)))) %>% `colnames<-`(colnames(log_or_vis_alt_df)) %>%
						mutate(Label = "Versus Secondary DENV2")
	ref_reg_label = as.data.frame(as.list(rep(NA, ncol(or_vis_alt_df))))%>% `colnames<-`(colnames(log_or_vis_alt_df))%>%
						mutate(Label = "Versus Asia")

	#Insert the reference seroprior label
	or_vis_alt_df = rbind(or_vis_alt_df[1, , drop = FALSE], ref_seroprior_label, or_vis_alt_df[-1, , drop = FALSE])
	log_or_vis_alt_df = rbind(log_or_vis_alt_df[1, , drop = FALSE], ref_seroprior_label, log_or_vis_alt_df[-1, , drop = FALSE])

	#Insert the reference reg label
	curr_n = nrow(or_vis_alt_df)
	or_vis_alt_df = rbind(or_vis_alt_df[1:(curr_n-1), , drop = FALSE], ref_reg_label, or_vis_alt_df[curr_n, , drop = FALSE])
	log_or_vis_alt_df = rbind(log_or_vis_alt_df[1:(curr_n-1), , drop = FALSE], ref_reg_label, log_or_vis_alt_df[curr_n, , drop = FALSE])

	# %%
	group_mean  = as.matrix(log_or_vis_alt_df[, c("PointEst_1997type", "PointEst_2009type", "PointEst_hospitalisation")])
	group_lower = as.matrix(log_or_vis_alt_df[, c("Lower_1997type", "Lower_2009type", "Lower_hospitalisation")])
	group_upper = as.matrix(log_or_vis_alt_df[, c("Upper_1997type", "Upper_2009type", "Upper_hospitalisation")])                   
	axis_lim = 5
	options(repr.plot.width = 12, repr.plot.height = 15)
	or_vis_plot = or_vis_alt_df %>% select(Label, DispVal_1997type, DispVal_2009type, DispVal_hospitalisation) %>% 
		forestplot(labeltext = c(Label, DispVal_1997type, DispVal_2009type, DispVal_hospitalisation),  
				mean = group_mean, lower = group_lower, upper = group_upper, clip = c(-axis_lim, axis_lim), 
				xticks = seq(-axis_lim, axis_lim, by = 1), ci.vertices = TRUE, boxsize = 0.1, xlab = "Est. Log Odds Ratio [95% CrI]",
				align = c("l", "l", "l", "r"), graphwidth = unit(2.5, "in"), zero = 0,
				legend = c("1997-type", "2009-type", "Hospitalisation")
				) |>
		fp_add_lines(h_2 = gpar(lty = 2)) |>
		fp_set_style(box = c("royalblue", "darkred", "darkgreen"), zero = "red", line = gpar(col = "black"),
				txt_gp = fpTxtGp(cex =1.2, legend = gpar(cex = 1.2), xlab = gpar(cex = 1.3), ticks = gpar(cex = 1))) |>
		fp_add_header(
					Label = "Variable",
					DispVal_1997type = "1997-type", 
					DispVal_2009type = "2009-type", 
					DispVal_hospitalisation = "Hospitalisation") |>
		
		fp_set_zebra_style("#EFEFEF")

	#save_plot(or_vis_plot, file.path(visuals_output_dir, paste0("OR_Group_Plot.png")), 12, 15, "in", 600)
	save_plot_all_formats(or_vis_plot, visuals_output_dir, "OR_Group_Plot", 12, 15, "in", 600)
}


# Versus Vaccine Trial Data ----
# Here, we compare the pooled probability of hospitalisation to the values from the TAK-003 clinical trials.
# %%
gen_vacc_comp_visual = function(merged_forest_df, visuals_output_dir, save_output = save_output){
	vacc_dir = file.path(base_dir, "Vaccine Trial Data")
	vacc_data_m57_sl = read_excel(file.path(vacc_dir, "Month57DataSriLanka.xlsx"))
	vacc_data_m57_sl = vacc_data_m57_sl %>% mutate(Scenario = paste0(Serotype, " ", Scope, " (Month 57)"), SET = "Month 57", TYPE = "DATA")
	vacc_data_m57_sl_bin = binconf(x = vacc_data_m57_sl$Hospitalised, n = vacc_data_m57_sl$VCD, method = "exact")
	vacc_data_m57_sl = vacc_data_m57_sl %>% cbind(vacc_data_m57_sl_bin) 
	vacc_data_m57_sl = vacc_data_m57_sl %>% select(-Scope) %>% rename(Cases = VCD, Severe = Hospitalised)

	pooled_ests = merged_forest_df %>% filter(SevClass == "hospitalisation") %>% #We are comparing trial data only to the pooled estimates where Prior Exposure is unknown
		separate(Scenario, sep = "-", into = c("Region", "PriorExposure", "Serotype"), remove = FALSE) %>%     
		filter(PriorExposure == "Unknown") %>% select(Serotype, N, Severe, Scenario, PointEst, Lower, Upper) %>% 
		rename(Cases = N) %>% mutate(SET = "POOLED", TYPE = "MODEL")

	vacc_comp = rbind(pooled_ests, vacc_data_m57_sl) %>% arrange(Serotype, SET)

	serotype_names = vacc_comp$Serotype %>% unique

	add_serotype_headers = function(curr_sero){
		curr_data = vacc_comp %>% filter(Serotype == curr_sero)
		curr_data = curr_data %>% mutate(PointEstDisp = round(PointEst * 100, 2), LowerDisp = round(Lower * 100, 2), UpperDisp = round(Upper * 100, 2)) %>% 
				mutate(DispText = sprintf("%.2f%% [%.2f%% to %.2f%%]", PointEst * 100, Lower * 100, Upper * 100)) %>% select(-c(PointEstDisp, LowerDisp, UpperDisp))
		curr_header = data.frame(Serotype = curr_sero, Cases = "", Severe = "", PointEst = NaN, Lower = NaN, Upper = NaN, Scenario = curr_sero, TYPE = "HEADER", SET = "HEADER", DispText = "")
		curr_data = curr_header %>% rbind(curr_data) %>% mutate(Scenario = ifelse(TYPE != "HEADER", paste0("     ", Scenario), Scenario))
		return(curr_data)
	}
	vacc_comp_vis_df = do.call(rbind, lapply(serotype_names, add_serotype_headers))
	vacc_comp_vis_df = vacc_comp_vis_df %>% mutate(Summ = SET == "POOLED") %>% 
	  mutate(DispText = str_replace_all(DispText, "\\.", "·"),
	         DispCases = ifelse(nchar(Cases) < 5, as.character(Cases), 
	                            str_replace(as.character(Cases), "(\\d)(?=(\\d{3})+$)", "\\1 "))) #Add a space in multiples of 3 digits
	
	vacc_comp_vis_df

	# %%
	options(repr.plot.width = 13, repr.plot.height = 14)
	vacc_comp_plot = vacc_comp_vis_df %>% 
	  forestplot(labeltext = c(Scenario, DispCases, Severe, DispText), is.summary = Summ, 
	             mean = PointEst, lower = Lower, upper = Upper,
	             xticks = seq(0, 1, by = 0.2), ci.vertices = TRUE, boxsize = 0.2,
	             align = c("l", "l", "l", "r", "r", "l"), graphwidth = unit(1.8, "in"), 
	  ) |>
	  fp_add_lines(h_2 = gpar(lty = 2)) |>
	  fp_set_style(box = "royalblue", line = gpar(col = "darkblue"), summary = "royalblue",
	               txt_gp = fpTxtGp(cex =1.2, ticks = gpar(cex = 1))) |>
	  fp_add_header(
	    Scenario = "Scenario",
	    DispCases = "N", 
	    Severe = "Hospitalised", 
	    DispText = "Proportion [95% CrI/CI]") |>
	  fp_set_zebra_style("#EFEFEF")
	
	#save_plot(vacc_comp_plot, file.path(visuals_output_dir, paste0("Vacc Placebo Comp.png")), 13, 14, "in", 600)
	save_plot_all_formats(vacc_comp_plot, visuals_output_dir, "Vacc Placebo Comp", 13, 14, "in", 600)
}

# Results Per Region and Country ----
# %%
generate_region_country_visual = function(model_inputs, index_df, visuals_output_dir, save_output = save_output){

	include_df = index_df %>% filter(FinalDecision %in% c("Include"))

	#If the CovidenceID is blank, use the CovidenceID_JanUpdate value (the JanUpdate IDs are new ones from the updated Covidence review we made)
	include_df = include_df %>% mutate(CovidenceID = ifelse(is.na(CovidenceID), CovidenceID_JanUpdate, CovidenceID))

	#Lookup table for the country that each study is in
	country_mapper = include_df %>% select(CovidenceID, Country) %>% distinct

	#Get the outcome data from the model input files
	outcomes_1997 = model_inputs[[1]]$outcome_df %>% mutate(SeveritySystem = "1997-type")
	outcomes_2009 = model_inputs[[2]]$outcome_df %>% mutate(SeveritySystem = "2009-type")
	outcomes_hosp = model_inputs[[3]]$outcome_df %>% mutate(SeveritySystem = "Hospitalisation")
	all_outcome_dfs = rbind(outcomes_1997, outcomes_2009, outcomes_hosp) %>% left_join(country_mapper, by = "CovidenceID")

	#Create a DataFrame that will serve as the basis of visualisations
	to_vis = all_outcome_dfs %>% select(CovidenceID, SeveritySystem, Region, Country) %>% distinct 

	#Find all studies with more than 1 country and split them into separate rows
	multi_country = to_vis %>% filter(str_detect(Country, ";")) %>% separate(Country, sep = ";", into = c("Country1", "Country2")) %>% 
		pivot_longer(c(Country1, Country2), values_to = "Country") %>% select(-name) %>% mutate(Country = trimws(Country))

	#Remove those with multiple countries and append the split row version
	to_vis = to_vis %>% filter(!str_detect(Country, ";")) %>% rbind(multi_country) %>% 
       		group_by(Region, Country, SeveritySystem) %>% 
			summarise(NumStudies = n_distinct(CovidenceID))

	options(repr.plot.width = 20, repr.plot.height = 16)
	to_vis %>% arrange(Country) %>% 
		ggplot(aes(x = Country, y = NumStudies, fill = SeveritySystem, label = NumStudies)) + 
		geom_bar(position = "stack", stat = "identity") +
		geom_text(size = 6, position = position_stack(vjust = 0.5), colour = "white") +
		facet_wrap(~Region, ncol = 1, scale = "free_x") +
		xlab("Location") + ylab("Number of Studies") +
		theme(text = element_text(size = 23),
				axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
				legend.justification = "left", legend.position = "top")
	ggsave(file.path(visuals_output_dir, "png", "StudyCountsSummary.png"), width = 20, height = 16, unit = "in")
	ggsave(file.path(visuals_output_dir, "eps", "StudyCountsSummary.eps"), width = 20, height = 16, unit = "in")
	ggsave(file.path(visuals_output_dir, "pdf", "StudyCountsSummary.pdf"), width = 20, height = 16, unit = "in")
}


#%%
visualise_model_set = function(model_set, effect_type, data_suffix = "", het_est = "stan", save_output = save_output){
	het_dir = file.path(base_dir, "Heterogeneity Estimates", io_set)
	visuals_output_dir = file.path(base_dir, "Visuals Output", io_set, model_set, "GroupVisuals")
  	print(visuals_output_dir)
	model_input_paths = file.path(model_input_dir, paste0("data_", sev_class_types, data_suffix, ".rds"))
	model_output_paths = file.path(model_output_dir, paste0("Results_", model_set, "_mean=0_sd=2_sd_mean=0.5_sdsd=2_", sev_class_types, ".rds"))

	# index_df = read_excel(file.path(index_dir, "ReviewIndex_Final.xlsx"), sheet = "Main")
	index_df = read_excel(file.path(index_dir, "ReviewIndex_Final_Strict.xlsx"), sheet = "Main")
	model_inputs = lapply(model_input_paths, readRDS)
	model_outputs = lapply(model_output_paths, readRDS)
	if(het_est == "separate"){
		het_results_paths = file.path(het_dir, paste0("I2_Estimates_", model_set, "_", sev_class_types, ".rds"))
		het_results = lapply(het_results_paths, readRDS)
	}else if(het_est == "stan"){
		het_results = NULL
	}else{
		stop("Invalid het_est value. Must be either 'stan' or 'separate'.")
	}
	
	#If visuals_output_dir does not exist, create it
	if(!dir.exists(visuals_output_dir)){
		dir.create(visuals_output_dir, recursive = TRUE)
	}


	merged_forest_df = gen_scenario_pooled_visual(model_inputs, model_outputs,sev_class_types, 
						het_results = het_results, visuals_output_dir = visuals_output_dir, effect_type = effect_type, save_output = save_output)

	gen_ref_or_table_visual(model_inputs, model_outputs, sev_class_types, visuals_output_dir = visuals_output_dir, save_output = save_output)
	gen_vacc_comp_visual(merged_forest_df, visuals_output_dir = visuals_output_dir, save_output = save_output)
	generate_region_country_visual(model_inputs, index_df, visuals_output_dir = visuals_output_dir, save_output = save_output)
}

#Function Calls ----
#%%
model_sets = c("LogisticRegression", "LogisticRegression_no_unknown", "FE", "NoCorr")
visualise_model_set("LogisticRegression", effect_type = "random", het_est = "stan", save_output = save_output)
visualise_model_set("LogisticRegression_no_unknown", effect_type = "random", data_suffix = "_no_unknown", het_est = "stan", save_output = save_output)
visualise_model_set("FE", effect_type = "fixed", het_est = "stan", save_output = save_output)
visualise_model_set("NoCorr", effect_type = "random", het_est = "stan", save_output = save_output)


#%%
# curr_input_data = model_inputs[[1]]
# curr_output_data = model_outputs[[1]]
# curr_sev_class_type = sev_class_types[[1]]
# curr_het_results = NULL
# #curr_input_data, curr_output_data, curr_sev_class_type, curr_het_results, effect_type = "random"
# model_set = "FE"
# effect_type = "fixed"
# het_est = "stan"



# het_dir = file.path(base_dir, "Heterogeneity Estimates", io_set)
# visuals_output_dir = file.path(base_dir, "Visuals Output", io_set, model_set, "GroupVisuals")
# model_input_paths = file.path(model_input_dir, paste0("data_", sev_class_types, ".rds"))
# model_output_paths = file.path(model_output_dir, paste0("Results_", model_set, "_mean=0_sd=2_sd_mean=0.5_sdsd=2_", sev_class_types, ".rds"))

# index_df = read_excel(file.path(index_dir, "ReviewIndex_Final.xlsx"), sheet = "Main")

# model_inputs = lapply(model_input_paths, readRDS)
# model_outputs = lapply(model_output_paths, readRDS)
# if(het_est == "separate"){
# 	het_results_paths = file.path(het_dir, paste0("I2_Estimates_", model_set, "_", sev_class_types, ".rds"))
# 	het_results = lapply(het_results_paths, readRDS)
# }else if(het_est == "stan"){
# 	het_results = NULL
# }else{
# 	stop("Invalid het_est value. Must be either 'stan' or 'separate'.")
# }


# #If visuals_output_dir does not exist, create it
# if(!dir.exists(visuals_output_dir)){
# 	dir.create(visuals_output_dir, recursive = TRUE)
# }

# #%%
# #Probabilities are given in the estimates of p 
# ref_region = curr_input_data$ref_region
# non_ref_region = ifelse(ref_region == "Asia", "Americas", "Asia")
# ref_serotype_exposure = curr_input_data$ref_serotype_exposure
# ref_scenario = curr_input_data$ref_scenario
# non_ref_seroprior = curr_input_data$char_mat_guide %>% filter(CharMatIndex >0) %>% pull(SeroPriorExp) %>% as.character


# scenario_labels = c(ref_scenario, 
# 				paste0(ref_region, "-", non_ref_seroprior),
# 				paste0(non_ref_region, "-", c(ref_serotype_exposure, non_ref_seroprior))
# 				)

# #Get estimates of p
# p_estimates = summary(curr_output_data, pars = "p")$summary %>% data.frame %>% 
# 			select(mean, X2.5., X97.5.) %>% rename(PointEst = mean, Lower = X2.5., Upper = X97.5.) %>% 
# 			mutate(Scenario = scenario_labels)

# scenario_df = curr_input_data$scenario_df %>% select(ScenIndex, Scenario, Inclusion, NumStudies, Severe, N)

# if(effect_type == "random"){
# 	tau_estimates = summary(curr_output_data, pars = "tau")$summary %>% data.frame %>% 
# 				select(mean, X2.5., X97.5.) %>% rename(PointEst = mean, Lower = X2.5., Upper = X97.5.) %>%
# 				mutate(ScenIndex = 1:nrow(.)) %>% mutate(DispTau = sprintf("%.2f [%.2f to %.2f]", PointEst, Lower, Upper))
				

# 	#If curr_het_results is NULL then use the I^2 estimates from Stan rather than the separate ones loaded from files
# 	if(is.null(curr_het_results)){
# 		curr_het_results = summary(curr_output_data, pars = "I2")$summary %>% data.frame %>% 
# 					select(mean, X2.5., X97.5.) %>% mutate(ScenIndex = 1:nrow(.)) %>% mutate(DispHet = sprintf("%.2f%% [%.2f to %.2f%%]", mean * 100, X2.5. * 100, X97.5. * 100))
# 		scenario_df = scenario_df %>% left_join(curr_het_results %>% select(ScenIndex, DispHet), by = "ScenIndex")
# 	}else{
# 		scenario_df = scenario_df %>% left_join(curr_het_results %>% select(Scenario, I2_est), by = "Scenario") %>%
# 					mutate(DispHet = ifelse(is.na(I2_est), "-", sprintf("%.2f%%", I2_est * 100))) %>% select(-I2_est)
# 	}
# 	#Join the tau estimates to the scenario_df then remvoe the ScenIndex column to simplify
# 	scenario_df = scenario_df %>% left_join(tau_estimates %>% select(ScenIndex, DispTau), by = "ScenIndex") 
# }

# #Get scenarios and then join the p estimates to them, removing rows where we do not have an estimate
# scenario_df = scenario_df %>% 
# 			left_join(p_estimates, by = "Scenario") %>% mutate(SevClass = curr_sev_class_type) %>%
# 			filter(!is.na(PointEst)) %>% arrange(Scenario) %>% 
# 			mutate(DispVal = sprintf("%.2f%% [%.2f to %.2f%%]", 100*PointEst, 100*Lower, 100*Upper)) %>% select(-ScenIndex) 

# disp_sev_class = ""
# if(curr_sev_class_type == "1997type"){
# 	disp_sev_class = "1997-type"
# }else if(curr_sev_class_type == "2009type"){
# 	disp_sev_class = "2009-type"
# }else{
# 	disp_sev_class = "Hospitalisation"
# }
# #Create the header row
# header_row = data.frame(Scenario = disp_sev_class, Inclusion = NA,
# 					NumStudies = NA, Severe = NA, N = NA, PointEst = NA, 
# 					Lower = NA, Upper = NA, DispVal = NA,
# 					SevClass = "HEADER")

# if (effect_type == "random"){
# 	header_row = header_row %>% mutate(DispHet = NA, DispTau = NA)
# }

# scenario_df = rbind(header_row, scenario_df %>% mutate(Scenario = paste0("       ", Scenario)))

# #%%
# if(is.null(het_results)){
# 	scenario_probs = mapply(process_scenario_probs, model_inputs, model_outputs, sev_class_types, 
# 						MoreArgs = list(curr_het_results = NULL, effect_type = effect_type), SIMPLIFY = FALSE)
# }else{
# 	#If het_results is not NULL, use those loaded from a file separately
# 	scenario_probs = mapply(process_scenario_probs, model_inputs, model_outputs, sev_class_types, het_results, 
# 						MoreArgs = list(effect_type = effect_type), SIMPLIFY = FALSE)
# }

# merged_scenarios = do.call(rbind, scenario_probs)
# merged_forest_df = merged_scenarios
# merged_forest_df = merged_forest_df %>% mutate(Pooled = ifelse(Inclusion == "Included", "Pooled", 
# 								ifelse(Inclusion == "Excluded", "Not Pooled", NA))) %>% 
# 									mutate(N = ifelse(Inclusion == "Excluded", "-", N),
# 									Severe = ifelse(Inclusion == "Excluded", "-", Severe))
