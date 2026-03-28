//Author: Hadrian Ang
//Date: 28/03/2026
/*
Overview
Code here is for a sensitivity analysis that uses a random effects model, but without the region-specific correlation matrices
for studies with multiple outcomes. This will let us test if region-specific differences are due to correlations or due to the 
region-specific fixed effects.
*/


data{
    int n_row; //Number of studies (also rows in total matrix)
    int n_col; //Number of serotype-prior exposures overall (usually 12, also cols in total matrix)
    int n_reg; //Number of regions (usually 2)
    int n_outcomes; //Number of outcomes, each a pair (n,e) 

    array[n_row, n_col] int total; //Number of dengue cases for the study in the i-th row with serotype-prior exposure j.
    array[n_row, n_col] int severe; //Number of severe cases for study in the i-th row with serotype-prior exposure j

    int n_sero_prior; //Number of serotype-prior exposures included in the analysis
    array[n_reg, n_col] int scenario_indices; //Look-up table to map region + serotype-prior exposure to an index (used in sigma)

    array[n_outcomes] int row_indices; //Row in the total matrix of the i-th outcome
    array[n_outcomes] int col_indices; //Col in the total matrix of the i-th outcome
    array[n_outcomes] int reg_indices; //Region index of the i-th outcome (usually, Americas: 1, Asia: 2)

    array[n_row] int region_vals; //Region index of the i-th study (usually, Americas: 1, Asia: 2)
    array[n_row] int num_results; //Number of outcomes contributed by the i-th study

    int num_scenarios; //Number of scenarios included in the analysis
    int reg_1_scen; //Number of included scenarios with region 1 (usually Americas)
    int reg_2_scen; //Number of included scenarois with region 2 (usually Asia)

    matrix[n_outcomes, (n_sero_prior + n_reg - 1)] char_matrix; //Characteristic matrix

    //Parameters for prior distributions

    //Prior for coefficients
    real coeff_prior_mean;
    real coeff_prior_sd; 

    //Prior for standard deviations for random effects
    real sd_prior_mean; 
    real sd_prior_sd;

    //Array to map scenario index to serotype-prior exposure index
    //We use this in computation of tau and tau^2
    array[num_scenarios] int scen_seroprior_ind_map;
    int n_tau_sim; //Number of simulations to estimate tau and tau^2
}

transformed data {
    //Array to map row and column to outcome index (which just goes from 1 to n_outcomes)
    array[n_row, n_col] int ind_mat; 
    for(i in 1:n_outcomes){
        ind_mat[row_indices[i], col_indices[i]] = i;
    } 

    //Number of pooled estimates we generate. Note that we can also generate pooled estimates for some excluded scenarios. 
    int n_reg_sero_prior = n_reg * n_sero_prior; 
    int char_mat_cols = n_sero_prior + n_reg - 1; //Number of columns in characteristic matrix

    //Characteristic matrix to generate summary at the end
    matrix[n_reg_sero_prior, char_mat_cols] summ_char_matrix = rep_matrix(0, n_reg_sero_prior, char_mat_cols);

    //Create a vector with the order of indices we want to create the summary char matrix
    array[n_reg_sero_prior] int summ_indices;
    summ_indices[1] = 0;
    summ_indices[n_sero_prior+1] = 0;
    for(i in 2:n_sero_prior){
        summ_indices[i] = i;
        summ_indices[i + n_sero_prior] = i;
    }
    
    //Characteristic matrix used to generate summary estimates for each scenario (including some excluded scenarios).
    //Each row is a scenario, each column is a characteristic (intercept, serotype-prior exposure, region).
    for(i in 1:rows(summ_char_matrix)){
        //Ones in the first column        
        summ_char_matrix[i, 1] = 1; 

        int sero_prior_ind = summ_indices[i]; 
        if(sero_prior_ind != 0){
            summ_char_matrix[i, sero_prior_ind] = 1; 
        }

        //If i is greater than 
        if(i > n_sero_prior){
            summ_char_matrix[i, char_mat_cols] = 1;
        }
    }

    int dims_beta_sero_prior = n_sero_prior - 1;
}

parameters{
    vector[dims_beta_sero_prior] beta_sero_prior; //Effect for serotype-prior exposure relative to reference
    real beta_reg; //Effect for region relative to reference
    real intercept; //The intercept, which also represents the reference scenario

    //Standard deviation of random-effects
    vector<lower = 0>[num_scenarios] sigma; 

    //Vector of random-effects. We draw a standard normal first, then transform to multivariate normal
    vector[n_outcomes] eps_base; 
}

transformed parameters {
    vector[n_outcomes] mu; //Fixed effects
    vector[n_outcomes] eps_vals; //Actual random-effect values
    vector[n_sero_prior + n_reg - 1] coeffs; //Coefficients for the logistic regression
    vector<lower = 0, upper = 1>[n_outcomes] theta; //Probabilities for the binomial

    //Stack of coefficients (starting with intercept, then serotype-prior exposures, and region as the last element)
    coeffs = append_row(intercept, append_row(beta_sero_prior, beta_reg));
    mu = char_matrix * coeffs; 

    //We transform the random-effects to their actual values
    //In this sensitivity analysis without correlations, we just multiply the base by the appropriate sigma.
    for(i in 1:n_outcomes){
        int k = reg_indices[i]; //Region index of outcome i
        int j = col_indices[i]; //Serotype-prior exposure index of outcome i
        eps_vals[i] = eps_base[i] * sigma[scenario_indices[k,j]]; //Transform base epsilon to actual epsilon by multiplying by the appropriate sigma value (based on region and serotype-prior exposure)
    }

    theta = inv_logit(mu + eps_vals); 
}

model{
    //Setting priors
    intercept ~ normal(coeff_prior_mean, coeff_prior_sd); 
    beta_sero_prior ~ normal(coeff_prior_mean, coeff_prior_sd); 
    beta_reg ~ normal(coeff_prior_mean, coeff_prior_sd); 
    sigma ~ normal(sd_prior_mean, sd_prior_sd); 
    eps_base ~ normal(0,1); 

    //Main likelihood computation / fitting
    for(i in 1:n_row){
        for(j in 1:n_col){
            if(total[i,j] > 0){
                severe[i,j] ~ binomial(total[i,j], theta[ind_mat[i,j]]);
            }
        }
    }
}

generated quantities {
    //We compute the main quantities of interest

    //The pooled estimates per scenario p 
    vector[n_reg_sero_prior] logit_p = summ_char_matrix * coeffs; 
    vector[n_reg_sero_prior] p = inv_logit(logit_p); 

    //Odds ratio form of the coefficients
    vector[dims_beta_sero_prior + 2] log_or_beta = coeffs; 
    vector[dims_beta_sero_prior + 2] or_beta = exp(log_or_beta); 

    //Pairwise odds ratios for serotype-prior exposure pairs
    matrix[dims_beta_sero_prior, dims_beta_sero_prior] log_sero_prior_or_ratios; 
    matrix[dims_beta_sero_prior, dims_beta_sero_prior] sero_prior_or_ratios; 
    for(i in 1:dims_beta_sero_prior){
        for(j in 1:dims_beta_sero_prior){
            log_sero_prior_or_ratios[i,j] = beta_sero_prior[i] - beta_sero_prior[j];
            sero_prior_or_ratios[i,j] = exp(log_sero_prior_or_ratios[i,j]); 
        }
    }
    
    //Pairwise odds ratio per scenario-pair
    matrix[n_reg_sero_prior, n_reg_sero_prior] log_scenario_or_ratios;
    matrix[n_reg_sero_prior, n_reg_sero_prior] scenario_or_ratios; 
    for(i in 1:n_reg_sero_prior){
        for(j in 1:n_reg_sero_prior){
            log_scenario_or_ratios[i,j] = logit_p[i] - logit_p[j]; 
            scenario_or_ratios[i,j] = exp(log_scenario_or_ratios[i,j]); 
        }
    }

    //We estimate tau and tau^2 - one for each scenario 
    vector[num_scenarios] tau2; 
    vector[num_scenarios] tau;

    //We estimate between-study variance through simulation for each included scenario.
    for(i in 1:num_scenarios){
        //Get the reg_sero_prior_ind for scenario i (mapping included scenarios to their row in summ_char_matrix / n_reg_sero_prior space)
        int reg_sero_prior_ind = scen_seroprior_ind_map[i]; 
        vector[n_tau_sim] eps_sim; //Simulated epsilon values for scenario i
        for(j in 1:n_tau_sim){
            eps_sim[j] = normal_rng(0, sigma[i]); //Simulate epsilon values based on the standard deviation of the random effects for scenario i
        }
        vector[n_tau_sim] theta_sim = inv_logit(logit_p[reg_sero_prior_ind] + eps_sim); //Simulate theta values based on the simulated epsilon values and the fixed effect for scenario i
        tau2[i] = variance(theta_sim); //Estimate tau^2 as the variance of the simulated theta values
        tau[i] = sd(theta_sim); //Get tau as the square root of tau^2
    }

    //We compute the prediction intervals for p, the pooled meta-analytic estimates for each scenario.
    //To get the prediction intervals, we just simulate a random effect value using sigma 
    vector[num_scenarios] pred_interval_p;
    for(i in 1:num_scenarios){
        pred_interval_p[i] = inv_logit(logit_p[scen_seroprior_ind_map[i]] + normal_rng(0, sigma[i]));
    }

}
