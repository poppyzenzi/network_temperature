library(dplyr)
library(tidyr)
library(boot)
library(psychonetrics)

## run remotely on Eddie cluster 
## bootstrap temperature CIs

# read in imputed data
abcd_mode_imputed <- read.table('abcd_mode_imputed.txt')
alspac_mode_imputed <- read.table('alspac_mode_imputed.txt')
mcs_mode_imputed <- read.table('mcs_mode_imputed.txt')

cohorts <- list(abcd = abcd_mode_imputed, alspac = alspac_mode_imputed, mcs = mcs_mode_imputed)

# symptom labels
labels_list <- list(
  alspac = c("unhappy", "anhedonia", "apathetic", "restless", "worthless","tearful", "distracted", "self_loathing", "guilty", "isolated", "unloved", "inadequate", "incompetent"),
  abcd = c("worthless","anxious","guilty","self_conscious","unhappy","worry"),
  mcs = c("malaise", "worries", "unhappy", "anxiety","fears")
)

# create empty list to store wide cohort dfs
wide_cohort_dfs <- list()

# iterate over cohorts > make wide dataframe > add to list
for (cohort_name in names(cohorts)) {
  cohort <- cohorts[[cohort_name]]
  wide_df <- reshape(as.data.frame(cohort) %>% dplyr::select(-sex), 
                     idvar = "id", timevar = "time", direction = "wide")
  wide_cohort_dfs[[cohort_name]] <- wide_df
}

# Function to bootstrap on individuals (wide) > convert to long > extract temperature
bootstrap_function <- function(data, indices, labels) {
  
  boot_data <- data[indices, ]
  boot_data$id <- seq(1, nrow(boot_data))
  
  boot_data_long <- boot_data %>%
    pivot_longer(
      cols = -id,
      names_to = c(".value", "time"),
      names_pattern = "(.*)\\.(\\d+)"
    )
  
  # Network analysis on boot data in long format
  model <- Ising(boot_data_long, vars = labels, groups = 'time') %>% groupequal("omega") %>% runmodel
  # Extract temperature estimates from model (1/beta)
  temp_est <- 1 / model@parameters$est[model@parameters$matrix == 'beta']
  # Omit t1 as it's fixed to 1 and stops boot.ci running
  temp_est <- temp_est[-1]
  
  return(temp_est)
}

# # bootstraps
num_bootstraps <- 1000

# empty list for bootstrap results
bootstrap_results <- list()

# loop to bootstrap cohorts
for (cohort_name in names(wide_cohort_dfs)) {
  wide_df <- wide_cohort_dfs[[cohort_name]]
  labels <- labels_list[[cohort_name]]
  bootstrap_results[[cohort_name]] <- boot(data = wide_df,
                                           statistic = bootstrap_function(data, indices, labels),
                                           R = num_bootstraps)
  saveRDS(bootstrap_results[[cohort_name]], paste0('bootstrap_results_', cohort_name, '.rds'))
}


# loop over cohorts and labels > perform bootstrap > save
for (cohort_name in names(wide_cohort_dfs)) {
  wide_df <- wide_cohort_dfs[[cohort_name]]
  labels <- labels_list[[cohort_name]]
  bootstrap_results[[cohort_name]] <- boot(
    data = wide_df,
    statistic = function(data, indices) bootstrap_function(data, indices, labels),
    R = num_bootstraps
  )
  saveRDS(bootstrap_results[[cohort_name]], paste0('bootstrap_results_', cohort_name, '.rds'))
}