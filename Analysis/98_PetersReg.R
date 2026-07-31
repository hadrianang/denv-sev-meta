# 97_PetersReg.R
#In this script, we just run Peters' regression for each of the three severity definitions.
#We use the implementation from the R meta package.
#%%
library(meta)
library(tidyverse)

#%%
io_set = "Main Results"
base_dir = getwd()
proc_data_dir = file.path(base_dir, "Processed Data", io_set)

sev_class_types = c("1997type", "2009type", "hospitalisation")
curr_sev_class = sev_class_types[3]
for(curr_sev_class in sev_class_types) {
  curr_input_name = paste0("data_", curr_sev_class, ".rds")

  curr_input_data = readRDS(file.path(proc_data_dir, curr_input_name))

  curr_outcome_df = curr_input_data$outcome_df

#Run Peters' regression for the current severity class type
  mp = metaprop(event = curr_outcome_df$Severe, n = curr_outcome_df$N, sm = "PLOGIT")
  bias_res = metabias(mp, method.bias = "peters")
  bias_df = data.frame(
    bias = bias_res$estimate["bias"],
    se.bias = bias_res$estimate["se.bias"],
    statistic = bias_res$statistic,
    df = bias_res$df,
    p.value = bias_res$p.value
  )
  rownames(bias_df) = NULL
  bias_df %>% write.csv(file.path(base_dir, "Visuals Output", io_set, "PublicationBias", paste0("PetersReg_", curr_sev_class, ".csv")), row.names = FALSE)
}
