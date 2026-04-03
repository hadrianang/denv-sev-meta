//Author: Hadrian Ang
//Date: 06/03/2025
/*
Overview
The code here is for a logistic regression model to be used for the meta-analysis
of dengue severity by serotype-prior exposure status and geographic region.

We define a scenario as some combination of region + serotype + prior exposure
*/

functions{
    //Function that takes a correlation matrix Rho and selects
    //relevant rows and columns based on an indicator vector
    //We use this to get only the part of the correlation matrix relevant to each
    //study i, so we can transform the random effects to a multivariate normal. 

    /*
    @params
    Rho: N x N correlation matrix corresponding to a given region r
    indic: N-d vector of numeric values. If indic[i] = 0, then we delete the i-th row and column in Rho. 
    num_vals: Integer denoting the number of non-zero elements in indic
    num_scens: Number of scenarios included in the analysis for the region r
    inclusion: N-d vector that we can use to see if the scenario corresponding 
                to some column in the totals matrix is included in the analysis. 
                If the element is less than 0, then the column is excluded in the analysis. 
    */
    matrix cut_submat(matrix Rho, vector indic, int num_vals, int num_scens, vector inclusion){
        //Selector matrix which we can multiply to Rho to get only the rows and columns we want
        matrix[num_vals, num_scens] selector = rep_matrix(0, num_vals, num_scens);
        
        int scenario_counter = 1;
        int counter = 1; 

        //For each element in the indicator vector
        for(i in 1:num_elements(indic)){
            //Note that if inclusion[i] < 0 then this scenario is not part of the analysis
            real curr_val = indic[i]; 
            real inc_check = inclusion[i]; 

            if(curr_val > 0){
                //If the value is positive, we are sure that the scenario is included
                //Set the selector element in the current row (counter) to 1
                selector[counter, scenario_counter] = 1;

                //Increment counter
                counter += 1; 
            }

            if(inc_check > 0){
                //If the scenario is included, we increment the scenario_counter
                scenario_counter += 1;
            }
        }
        matrix[num_vals, num_vals] to_ret = selector * Rho * selector'; 

        return to_ret; 
    }
    //Alternative version of taking appropriate sub-matrix of Rho using the same parameters
    matrix cut_submat_alt(matrix Rho, vector indic, int num_vals, int num_scens, vector inclusion){
        //Matrix to return
        matrix[num_vals, num_vals] to_ret;
        array[num_vals] int to_get; //row and column indices of Rho that should be returned
        
        int curr_get_ind = 1; //Index in the to_get array we are on
        int rho_ind = 1; //Index in Rho we are on (Rho does not match the dimensions in indic since some scenarios are excluded)
        for(i in 1:num_elements(indic)){
            //If indic[i] > 0 then the returned matrix should contain the row and column of Rho that we are currently on (rho_ind)
            if(indic[i] > 0){
                to_get[curr_get_ind] = rho_ind; //Add rho_ind to the selector array
                curr_get_ind += 1; //Increment index on the to_get array for when we place the next value
            }
            if(inclusion[i] > 0){
                rho_ind += 1; //If the scenario that we just passed is included in the analysis, increment the index on the Rho matrix
            }
        }

        to_ret = Rho[to_get, to_get]; //Use the index array to retrieve the required rows and columns of Rho

        return to_ret; 
    }
}


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

    //Array to map scenario index to reg + serotype-prior exposure index
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

    //Correlation matrix for each region (size depends on number of included scenarios in each)
    corr_matrix[reg_1_scen] Rho_1; 
    corr_matrix[reg_2_scen] Rho_2; 
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

    //For each study in the analysis
    for(i in 1:n_row){
        int k = region_vals[i]; //Region index of study i
        int num_res = num_results[i]; //Number of outcomes from study i

        //If study i gives a single outcome,
        //then we convert the random effect value to eps ~ N(0, sigma) by multiplying sigma the base standard normal value
        if(num_res == 1){
            for(j in 1:n_col){
                if(total[i,j] > 0){
                    eps_vals[ind_mat[i,j]] = eps_base[ind_mat[i,j]] * sigma[scenario_indices[k,j]];
                    //print("eps_vals SINGLE: ",  eps_vals[ind_mat[i,j]], ", OUTCOME: ", ind_mat[i,j]); 
                    break; 
                }
            }
        }else{
            //If study i has multiple outcomes, then we convert to eps ~ MvNormal(0, Sigma)
            vector[num_res] eps_temp; //eps_base values relevant to study i
            vector[num_res] sigma_vec; //sigma values relevant to study i

            matrix[num_res, num_res] Sub_rho; //Sub-matrix of the correlation matrix Rho containing rows and columns relevant to study i
            
            int start_ind = -1; //First index in eps_vals relevant to study i (outcome space)
            int counter = 1; //Collector index in num_res space
            for(j in 1:n_col){
                //If there is a valid outcome (from study i)
                if(total[i,j] > 0){
                    //Collect the vector of base epsilon values
                    eps_temp[counter] = eps_base[ind_mat[i,j]]; //Get the epsilon value corresponding to the outcome i,j
                    if(start_ind == -1){
                        start_ind = ind_mat[i,j]; //Start of the epsilon indices relevant to study i
                    }
                    //Collect relevant standard deviations
                    sigma_vec[counter] = sigma[scenario_indices[k,j]]; //Get the value of sigma for the current scenario (based on reigon and serotype-prior exposure)
                    counter += 1; //Increment collector index

                    //Stop if we already have all the results from i
                    if(counter == (num_res+1)){
                        break;
                    }
                }
            }

            //Build correlation matrix
            //If study i is in region 1 (usually Americas)
            if(k==1){
                Sub_rho = cut_submat_alt(Rho_1, to_vector(total[i,]), num_res, reg_1_scen, to_vector(scenario_indices[k,]));
            }else{
                //Assumes only 2 regions, thus these are for studies in region 2 (usually Asia)
                Sub_rho = cut_submat_alt(Rho_2, to_vector(total[i,]), num_res, reg_2_scen, to_vector(scenario_indices[k,]));
            }
            vector[num_res] eps_trans = sigma_vec .* (cholesky_decompose(Sub_rho) * eps_temp); //transform the base eps values to the multivariate normal
            eps_vals[start_ind: (start_ind + num_res - 1)] = eps_trans; //set the vector storing the transformed eps values
            //print("eps_temp: ", eps_temp);
            //print("eps_trans: ", eps_trans, ", OUTCOMES: ", start_ind, " TO ", (start_ind + num_res - 1));
        }
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
    Rho_1 ~ lkj_corr(2);
    Rho_2 ~ lkj_corr(2);

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

    //We can give an estimate of I^2 by first computing v_x, the binomial variance of outcome x 
    vector[n_outcomes] v_x;
    vector[num_scenarios] v_sums = rep_vector(0, num_scenarios); //This is the variance of the pooled estimate for the scenario that outcome x belongs to, which we can get from the generated p values and the total number of cases for outcome x. We use this as an approximation for the between-study variance since we don't have a closed form for the variance of the predicted probabilities.
    vector[num_scenarios] v_counts = rep_vector(0, num_scenarios);
    for(i in 1:n_outcomes){
        int curr_scenario = scenario_indices[reg_indices[i], col_indices[i]]; //Scenario index of outcome i
        real curr_p = p[scen_seroprior_ind_map[curr_scenario]]; // Get the value of the pooled estimate for the scenario that outcome i belongs to

        v_x[i] = (curr_p * (1 - curr_p)) / total[row_indices[i], col_indices[i]]; //Binomial variance for outcome i
        v_sums[curr_scenario] += v_x[i]; //Add the variance of outcome i to the sum of variances for the scenario that outcome i belongs to
        v_counts[curr_scenario] += 1; //Increment the count of outcomes for the scenario that outcome i belongs to
    }

    vector[num_scenarios] v_bar = v_sums ./ v_counts; //Get the average variance for each scenario
    vector[num_scenarios] I2 = tau2 ./ (tau2 + v_bar); //Estimate I^2 using the formula I^2 = tau^2 / (tau^2 + v_bar)
    


}
