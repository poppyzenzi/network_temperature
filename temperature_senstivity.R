################################################################################
########################## SENSITIVITY ANALYSES ################################

## This script gives examples to derive network entropy, connectivity and symptom variance
## Downstream analyses for the calculation of network temperature
## We continue the example with the ALSPAC cohort

# load required packages
library(dplyr)
library(ggplot2)
library(IsingSampler)
library(NetworkComparisonTest)

# load data
best.model.alspac <- ('/path/to/saved/model/output')
alspac_mode_imputed <- ('/path/to/imputed/symptom/data')


################## 1. Calculate network entropy #####################

# create list to store values and specify cohort timepoints and ages
entropy_values <- list()
time_points <- c(1:6)
ages <- c(11,13,14,17,18,19)

# extract edge weights (omega), thresholds (tau), and temp (beta) from each model
for (t in time_points) {
  graph <- best.model.alspac@modelmatrices[[t]]$omega
  thresholds <- best.model.alspac@modelmatrices[[t]]$tau
  beta <- as.vector(best.model.alspac@modelmatrices[[t]]$beta)
  
  if (!isSymmetric(graph)) {
    graph <- (graph + t(graph)) / 2  # Ensure symmetry if necessary
  }
  
  # store in list of entropy values
  entropy <- IsingSampler::IsingEntrophy(graph = graph, thresholds = thresholds, beta = beta)
  entropy_values <- c(entropy_values, entropy)
}

# set up plot data
plot_entropy_data <- data.frame(Age=ages,
                                Entropy=unlist(entropy_list))

# plot entropy change with age
ggplot(plot_entropy_data, aes(x = Age, y = Entropy)) +
  geom_point(color= "blue") +
  geom_line(color = "blue") +
  theme_minimal() +
  scale_x_continuous(breaks = ages) +
  labs(title = "ALSPAC", x = "Age", y = "Entropy") +
  theme(legend.position = "bottom")



############ 2. Derive/compare with global connectivity change ############

# Using the NetworkComparisonTest package
# recode the imputed data to 0,1 (required by NCT)
alspac_mode_imputed_recoded <- alspac_mode_imputed %>%
  mutate(across(4:ncol(alspac_mode_imputed), ~case_when(
    . == 1 ~ 1,
    . == -1 ~ 0,
    TRUE ~ .
  )))

# Extract data for each time point
time_data <- list()
for (t in 1:6) {
  time_data[[t]] <- alspac_mode_imputed_recoded %>%
    filter(time == t) %>%
    dplyr::select(all_of(vars)) %>%
    filter(!if_all(everything(), is.na))
}

# vector to store strength estimates
strength_estimates <- numeric(length = 6)

# Perform Network Comparison Test between consecutive time points
for (t in 1:(length(time_data) - 1)) {
  nct_result <- NCT(
    data1 = time_data[[t]],
    data2 = time_data[[t + 1]],
    it = 200, # Use more iterations for real analysis (e.g., 1000+, do remotely on HPC for memory)
    binary.data = TRUE
  )
  
  # Save the strength estimates for each time point
  if (t == 1) { 
    # For the first pair, save both time points' strength estimates
    strength_estimates[t] <- nct_result$glstrinv.sep[1]
    strength_estimates[t + 1] <- nct_result$glstrinv.sep[2] 
  } else {
    # For subsequent pairs, save only the second time point's strength estimate
    strength_estimates[t + 1] <- nct_result$glstrinv.sep[2]
  }
}

# Set up plotting data
strength_estimates_df <- data.frame(Age = ages,StrengthEstimate = strength_estimates)

# Plot the strength estimates across age
ggplot(strength_estimates_df, aes(x = Age, y = StrengthEstimate)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Global Strength Estimate Over Time",
    x = "Time Point",
    y = "Global Strength Estimate"
  ) +   ylim(0, 50) + scale_x_continuous(breaks = ages) +
  theme_minimal()


############# 3. Calculate symptom mean and variance across time ##############

# function to calculate mean and variances of each column (symptom)
calculate_mean <- function(matrix_data) {
  apply(matrix_data, 2, function(column) {
    mean(column, na.rm = TRUE) 
  })
}

calculate_variance <- function(matrix_data) {
  apply(matrix_data, 2, function(column) {
    var(column, na.rm = TRUE)
  })
}

# Define the time points and lists
time_points <- c(1:6)
ages <- c(11,13,14,17,18,19)
vars <- names(alspac_mode_imputed)[4:ncol(alspac_mode_imputed)]
variances <- list() # list for symptom variances
mean_variances <- numeric() # list for average symptom variances across time
means <- list() # list for symptom means
overall_means <- numeric() # list for overall symptom means across time

# Loop through the time points and calculate variances
for (t in time_points) {
  time_data <- alspac_mode_imputed %>%
    filter(time == t) %>%
    dplyr::select(4:ncol(alspac_mode_imputed)) %>%
    as.matrix()
  
  variances[[as.character(t)]] <- calculate_variance(time_data)
  mean_variances <- c(mean_variances, mean(variances[[as.character(t)]]))
  
  means[[as.character(t)]] <- calculate_mean(time_data)
  overall_means <- c(overall_means, mean(means[[as.character(t)]]))
}

# Create a comparison data frame for the mean variances
comparison <- data.frame(
  Matrix = paste("Matrix", time_points, sep=""),
  MeanVariance = mean_variances,
  OverallMean = overall_means
)

print(comparison)

# Option to plot results overall symptoms
plot_all_data <- data.frame(Variance = mean_variances, Mean = overall_means, Age = factor(ages))
ggplot(plot_all_data, aes(x = Age, y = Mean)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = Mean-Variance, ymax = Mean+Variance), width = 0.1) +
  theme_minimal() +
  labs(title = "ALSPAC", x = "Age", y = "Mean symptom score") +
  theme(legend.position = "bottom")
