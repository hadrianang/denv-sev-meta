//Author: Hadrian Ang
//Date: 25/04/2025
/*
Overview
Code here is a simplified logistic regression that excludes the correlations in random effects 
and region-specific effects. In this case, we are mainly only fitting the serotype-prior exposure effects, random-effects,
and the variance of random effects.
*/

data{
    int n_row;
    int n_col;
    int n_outcomes;

    array[n_row, n_col] int total;
    array[n_row, n_col] int severe;

    array[n_outcomes] int row_indices;
    array[n_outcomes] int col_indices;
    
    array[n_row] int num_results;

    int n_sero_prior;

    matrix[n_outcomes, n_sero_prior] char_matrix;

    //Parameters for prior distributions

    //Prior for coefficients
    real coeff_prior_mean;
    real coeff_prior_sd; 

    //Prior for standard deviations for random effects
    real sd_prior_mean; 
    real sd_prior_sd; 
}

transformed data {
    array[n_row, n_col] int ind_mat;
    for(i in 1:n_outcomes){
        ind_mat[row_indices[i], col_indices[i]] = i;
    }

    matrix[n_sero_prior, n_sero_prior] summ_char_matrix;
    summ_char_matrix = diag_matrix(rep_vector(1, n_sero_prior));
    summ_char_matrix[,1] = rep_vector(1, n_sero_prior);
    int dims_beta_sero_prior = n_sero_prior - 1; 
}

parameters{
    vector[dims_beta_sero_prior] beta_sero_prior;
    real intercept;
}

transformed parameters {
    vector[n_outcomes] mu;
    vector[n_sero_prior] coeffs;
    vector<lower = 0, upper = 1>[n_outcomes] theta; 

    coeffs = append_row(intercept, beta_sero_prior); 
    mu = char_matrix * coeffs; 

    theta = inv_logit(mu);
}

model{
    //Setting priors
    intercept ~ normal(coeff_prior_mean, coeff_prior_sd); 
    beta_sero_prior ~ normal(coeff_prior_mean, coeff_prior_sd); 
    
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
    //The pooled estimates per scenario p 
    vector[n_sero_prior] logit_p = summ_char_matrix * coeffs; 
    vector[n_sero_prior] p = inv_logit(logit_p); 
    vector[n_sero_prior] log_or_beta = coeffs; 
    vector[n_sero_prior] or_beta = exp(log_or_beta); 

    //OR ratios
    matrix[dims_beta_sero_prior, dims_beta_sero_prior] log_sero_prior_or_ratios; 
    matrix[dims_beta_sero_prior, dims_beta_sero_prior] sero_prior_or_ratios; 
    for(i in 1:dims_beta_sero_prior){
        for(j in 1:dims_beta_sero_prior){
            log_sero_prior_or_ratios[i,j] = beta_sero_prior[i] - beta_sero_prior[j]; 
            sero_prior_or_ratios[i,j] = exp(log_sero_prior_or_ratios[i,j]); 
        }
    }
}
