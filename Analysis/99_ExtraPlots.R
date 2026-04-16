#99_ExtraPlots.R
#The code here is primarily for generating some extra plots that we can use for our peer review response. 

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

#%%
curr_sev_class_type = sev_class_types[1]
