# Association of Dengue Severity with the Serotype of Infection: a systematic review and meta-analysis
### Authors  
Hadrian Jules Ang<sup>1</sup>   
Clare McCormack<sup>1</sup>  
Christl A. Donnelly<sup>1,2</sup>  
Ilaria Dorigatti<sup>1</sup>  

<sup>1</sup>Medical Research Council Centre for Global Infectious Disease Analysis, School of Public Health, Imperial College London, London, United Kingdom

<sup>2</sup>Department of Statistics, University of Oxford, Oxford, United Kingdom

## Directory Structure
<pre>
denv-sev-meta/
├── Analysis/                                  
│   ├── Cohort Analysis/                       #Analysis for cohort study sub-analysis    
│   │   ├── Model Inputs/                      #Input for sub-analysis model     
│   │   ├── Model Outputs/                     #Output of sub-analysis model    
│   │   ├── Output Visuals/                    #Plots of results for sub-analysis   
│   │   ├── 6_Data Processing.ipynb            #Prepares data for sub-analysis model     
│   │   ├── 7_RunCohortModel.R                 #Runs the sub-analysis model
│   │   ├── 8_CohortVisuals.ipynb              #Generates plots of sub-analysis results
│   │   └── CohortLogisticRegression_FE.stan   #Defines sub-analysis Stan model  
│   ├── Heterogeneity Estimates/               #Contains heterogeneity estimates  
│   ├── Model Code/                            #Stan code for main analysis model  
│   ├── Model Output/                          #Results of main model for 3 severity systems
│   ├── Processed Data/                        #Processed data for input into main model
│   ├── Vaccine Trial Data/                    #Contains vaccine trial data for validation 
│   ├── Visuals Output/                        #Output plots from main analysis
│   ├── 1_Data Preparation.ipynb               #Code prepares data for input into model
│   ├── 2_RunModel.R                           #Runs the main model  
│   ├── 3_Heterogeneity Computations.ipynb     #Estimates heterogeneity using model results - old version, now computed in Stan
│   ├── 4_Severity Class Visuals.R             #Generates visuals for a single severity system
│   ├── 5_Group Visuals.R                      #Generates some joint visuals across systems
│   ├── 6_ComparisonPlots.R                    #Generates some plots to compare results across model runs
│   ├── 7_DiagPlots.R                          #Generates trace plots and tables of parameter estimates
│   ├── 98_PetersReg.R                         #Run standard Peters' regression for assessing small study effects
│   └── 99_ExtraPlots.R                        #Generates additional plots such as those for assessing small study effects
├── Data/                                      #Extracted data for meta-anaylsis  
├── ReviewIndex_Final_Strict.xlsx              #Index to link study information to extracted data - for sensitivity analysis
├── ReviewIndex_Final.xlsx                     #Index to link study information to extracted data
└── README.md  
</pre>

## Packages and Versions Used

R packages:
1. tidyverse = 2.0.0
2. ggplot2 = 4.0.2
3. rstan = 2.32.7
4. Hmisc = 5.2-5
5. readxl = 1.4.5
6. writexl = 1.5.4

In addition, we used Stan version 2.32.2.

## Usage
With the exception of the main and sub-analysis model coded in Stan, the rest of the code in this repository is written in R. The `ReviewIndex_Final.xlsx` file contains some information extracted from each study (location, date, age of participants, serotyping method, prior exposure statuses found, severity guidelines applied), while also linking each study to an Excel file containing extracted data on the number of cases with each severity classification per serotype, per prior exposure status found (kept in `Data/` folder). 

Code should be run in the order given by the prefix number. For example, `1_Data Preparation.ipynb` should be run before `2_RunModel.R`. 
1. `1_Data Preparation.ipynb` - reads the index and extracted case data and outputs RDS files containing input necessary for the Stan analysis models into the `Processed Data/` folder.
2. `2_RunModel.R` - contains code for actual Stan sampling using the model defined in the `Model Code` folder for one severity classification system. The resulting stanfit object is then written to the `Model Output/` folder. This has to be run three times, once each for 1997-type, 2009-type, and hospitalisation severity. 
3. `3_Heterogeneity Computations.ipynb` - computes estimates of heterogeneity per scenario as measured by $I^2$. Results are written to the `Heterogeneity Estimates/` folder and this should again be run three times, once for each severity classification system. However, results from this are an old version, and the current heterogeneity estimates are computed directly in Stan. 
4. `4_Severity Class Visuals.ipynb` - code generates visuals for a single severity classification system written to the `Visuals Output/Main Results/<SeveritySystem>/` folder. For example to `Visuals Output/Main Results/1997type/` for 1997-type severity. This code has to again be run once for each of the three severity classification systems.
5. `5_Group Visuals.R` - code generates visuals that involve all three severity classification systems and outputs them to the `Visuals Output/Main Results/GroupVisuals/` folder.
6. `6_ComparisonPlots.R` - reads results of model runs from `Model Outputs` then generates comparison plots in `Visuals Output/Main Results/ComparisonPlots/`
7. `7_DiagPlots.R` - generates trace plots and parameter estimate tables in `Visuals Output/Main Results/<SeveritySystem>/DiagPlots/`
8. `98_PetersReg.R` - code for running standard Peters' regression using the meta package. 
9. `99_ExtraPlots.R` - generates additional visuals, primarily for assessing small study effects.
10. `6_Data Processing.ipynb` - pre-processes data for input into the cohort sub-analysis model and outputs the results to the `Cohort Analysis/Model Inputs/` folder.
11. `7_RunCohortModel.R` - contains code for Stan sampling using the cohort sub-analysis model defined in `CohortLogisticRegression_FE.stan`. The resulting stanfit object is written to the `Cohort Analysis/Model Outputs/` folder.
12. `8_CohortVisuals.ipynb` - generates visuals for the results of the cohort study sub-analysis and writes them to the `Cohort Analysis/Output Visuals/` folder. 