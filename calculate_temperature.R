############################# Network Temperature ##############################

## This is a script to derive network temperature of depression symptoms using multi-group Ising models where group = time
## This is an example using the ALSPAC cohort

# load required packages 
library(bootnet)
library(psychonetrics)
library(dplyr)
library(qgraph)
library(mice)
library(ggmice)
library(VIM)
library(boot)
library(purrr)

############################## 1. DATA AND QC ##################################
# check/set your working directory first
getwd()

# read in symptom level data
smfq_symptoms <- read.table('smfq_symptoms_qcd.txt', check.names = FALSE)

# add any environmental variables required (e.g. for imputation or validation)
env <- read.table('alspac_envi_vars.txt', check.names=FALSE)
smfq_qcd <- merge(smfq_symptoms, env, by='id') %>% rename(`sex` = kz021)

# variable labels
labels <- c("unhappy", "anhedonia", "apathetic", "restless", "worthless",
            "tearful", "distracted", "self_loathing", "guilty", "isolated", 
            "unloved", "inadequate", "incompetent")
env_labels = list(names(smfq_qcd[,16:ncol(smfq_qcd)]))
colnames(smfq_qcd) <- c('id', 'time', unlist(labels), unlist(env_labels))

# recode variables [-1,1]
smfq_qcd <- smfq_qcd %>%
  mutate(across(3:15, ~case_when(
    . == 1 ~ 1,
    . == 0 ~ -1,
    TRUE ~ .
  )))

############################## 2. IMPUTATION ##################################
# Prepare the df to impute into and merge with existing data
template <- expand.grid(id = unique(smfq_qcd$id), time = 1:6)
prepped_df <- merge(template, smfq_qcd, by = c("id", "time"), all.x = TRUE)

# sort df by id and time
prepped_df <- prepped_df %>% arrange(id, time)

# ensure NA for non-time invariant columns
prepped_df <- prepped_df %>%
  mutate(across(unhappy:incompetent, ~ifelse(. %in% c("", "NA"), NA, .)))

# fill sex (have checked that individuals don't have conflicitng sex)
prepped_df <- prepped_df %>%
  group_by(id) %>%
  fill(sex, .direction = "down") %>% fill(sex, .direction = "up")

# reorder cols
prepped_df <- prepped_df[, c("id", "time", "sex", unlist(labels))]

# Check completed df
print(prepped_df)

# inspect missing data patterns (2 options)
plot_pattern(prepped_df, square = TRUE, rotate = TRUE) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("ALSPAC")
aggr_plot <- aggr(prepped_df, col=c('navyblue','red'), numbers=TRUE, 
                  sortVars=TRUE,labels=names(prepped_df), cex.axis=.7, gap=3, 
                  ylab=c("Histogram of missing data","Pattern"))

# make predictor matrix for symptom vars
predMat <- make.predictorMatrix(prepped_df)
predMat[,1] <- 0
meth <- make.method(prepped_df)
meth[c("sex")] <- ""

# impute (takes a while, you might want to do this remotely on an HPC)
imputed_smfq <- mice(prepped_df,m=41,maxit=20,meth=meth,seed=123,
                     predictorMatrix=predMat)
saveRDS(imputed_smfq, 'imputed_smfq_41.rds')

# extract and combine pooled datasets into long format
imputed_list_smfq <- lapply(1:41, function(i) complete(imputed_smfq, action = i))
all_imputed_smfq <- bind_rows(imputed_list_smfq)

# sanity check: count the number of duplicate rows based on 'id' and 'time' (here, 9217*6 = 55302)
nrow(all_imputed_smfq %>%
       group_by(id, time, sex) %>%
       filter(n() > 1) %>%
       summarise(count = n()) %>%
       ungroup())

# group by 'id', 'time', 'sex' and the symptom columns, then count occurrences
combination_counts <- all_imputed_smfq %>%
  group_by(id,time,sex,unhappy,anhedonia,apathetic,restless,worthless,tearful,
           distracted,self_loathing,guilty,isolated,unloved,inadequate,incompetent) %>%
  summarise(count = n(), .groups = 'drop')

# for each grouping of 'id', 'time', 'sex', select the combination with the highest count
alspac_mode_imputed <- combination_counts %>%
  group_by(id, time, sex) %>%
  slice_max(count, n = 1, with_ties = FALSE) %>%
  ungroup() %>% dplyr::select(-count)

# final dataset that pools the mode symptom combination from multiple imputation
head(alspac_mode_imputed)
# save
write.table(alspac_mode_imputed, 'alspac_imputed_ising_mode.txt', col.names=TRUE)


########################## 3. MULTIGROUP ISING MODEL ############################

# list the nodes in the network
vars <- names(alspac_mode_imputed)[4:ncol(alspac_mode_imputed)]

# Form saturated model and run [all params free] (estimator = 'ML')
model1 <- Ising(alspac_mode_imputed, vars = vars, groups = "time") %>% runmodel
# Prune-stepup to find a sparse model:
model1b <- model1 %>% prune(alpha = 0.05) %>%  stepup(alpha = 0.05)
# Equal networks (omega = network structure, edges btw nodes equal across time):
model2 <- model1 %>% groupequal("omega") %>% runmodel
# Prune-stepup to find a sparse model:
model2b <- model2 %>% prune(alpha = 0.05) %>% stepup(mi = "mi_equal", alpha = 0.05)
# Equal thresholds (tau = threshold/intercept structure, omega still equal plus external fields equal):
model3 <- model2 %>% groupequal("tau") %>% runmodel
# Prune-stepup to find a sparse model:
model3b <- model3 %>% prune(alpha = 0.05) %>% stepup(mi = "mi_equal", alpha = 0.05)
# Equal beta (beta = inverse temperature, omega and tau also constrained):
model4 <- model3 %>% groupequal("beta") %>% runmodel
# Prune-stepup to find a sparse model:
model4b <- model4 %>% prune(alpha = 0.05) %>% stepup(mi = "mi_equal", alpha = 0.05)

# Compare all models, increasing constraints (RMSEA <0.05 good, minimise AIC and BIC)
psychonetrics::compare(
  `1. all parameters free (dense)` = model1, # saturated
  `2. all parameters free (sparse)` = model1b,
  `3. equal networks (dense)` = model2, # equal network structure 
  `4. equal networks (sparse)` = model2b,
  `5. equal networks and thresholds (dense)` = model3, # structure and external fields
  `6. equal networks and thresholds (sparse)` = model3b,
  `7. all parameters equal (dense)` = model4, # structure, external and temperature 
  `8. all parameters equal (sparse)` = model4b
) %>% arrange(BIC) 

# select best model manually based on combination of model parameters and wider context
# here the best model is model2
best.model.alspac <- model2

# extract edge weights and plot the network
all_network_smfq <- getmatrix(best.model.alspac, "omega")
network_smfq <- getmatrix(best.model.alspac, "omega")[[1]]
graph_smfq <- qgraph(network_smfq, layout = 'spring', labels = vars, 
                     theme = 'colorblind', label.prop=0.99, node.width=1.4, 
                     label.norm='000000')

# (optional) centrality plot
centralityPlot(list(Wave1 = all_network_smfq[[1]]),
               theme_bw=FALSE, scale = "z-scores", 
               include = c("Strength","Closeness","Betweenness"), 
               labels = vars, orderBy = 'Strength')


# another nice visualisation: plot heatmaps of symptom covar matrix
rownames(network_smfq) <- labels[1:13]
colnames(network_smfq) <- labels[1:13]

heatmap(network_smfq, 
        symm = TRUE,
        col = viridis::plasma(100),
        Rowv = NA,
        main = paste("SMFQ items"))

# extract temperature and calculate 95%CIs from standard errors of the beta parameter
temp_smfq <- as.numeric(lapply(getmatrix(best.model.alspac, "beta"), 'mean'))
betas_est <- best.model.alspac@parameters$est[best.model.alspac@parameters$matrix == "beta"]
betas_se <- best.model.alspac@parameters$se[best.model.alspac@parameters$matrix == "beta"]
z = qnorm(0.975)
upperCI <- temp_smfq + (z*betas_se)
lowerCI <- temp_smfq - (z*betas_se)

# set up the plotting data
ages <- c(11,13,14,17,18,19)
alspac_temp_data <- data.frame(Age = ages,
                                     Temperature = 1/temp_smfq,
                                     UpperCI = 1/upperCI,
                                     LowerCI = 1/lowerCI)


# plot temperature change across age
ggplot(alspac_temp_data, aes(x = Age, y = Temperature)) +
  geom_point(shape = 16) +
  geom_line() +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.1) +
  scale_x_continuous(breaks = ages) +
  labs(x = "Age", y = "Temperature", title = "Temperature change") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 12), legend.position='none')


####################### BOOTSTRAPPING #############################

# read in bootstrap results from HPC
# see 'boot.ci' script
boot.alspac <- readRDS('bootstrap_results_alspac.rds')
head(boot.alspac)

# get bootstrapped CI's
std_errors <- c(NA,apply(boot.alspac$t, 2, sd))
z = qnorm(0.975)
upper_boot_CI <- temp_smfq + (z*std_errors)
lower_boot_CI <- temp_smfq - (z*std_errors)

# set up plot data
boot_ci_data <- data.frame(Age = ages,
                           Temperature = 1/temp_smfq,
                           UpperCI = 1/upperCI,
                           LowerCI = 1/lowerCI,
                           UpperBootCI = 1/upper_boot_CI,
                           LowerBootCI = 1/lower_boot_CI)

boot_ci_data <- boot_ci_data %>%
  mutate(CI_Type = 'Analytical 95% CIs',
         BootCI_Type = 'Bootstrapped 95% CIs')

# plot bootstrapped CIs
ggplot(boot_ci_data, aes(x = Age, y = Temperature)) +
  geom_point(shape = 16) +
  geom_line() +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI, color = CI_Type), width = 0.1, show.legend = TRUE) +
  geom_errorbar(aes(ymin = LowerBootCI, ymax = UpperBootCI, color = BootCI_Type), width = 0.1, show.legend = TRUE) +
  labs(x = "Age (years)", y = "Temperature") +
  theme_minimal() +
  theme(axis.ticks.length = unit(14, "points")) +
  scale_x_continuous(breaks = ages) +
  scale_color_manual(values = c('Analytical 95% CIs' = 'purple', 'Bootstrapped 95% CIs' = 'red'), name=NULL) +
  guides(color = guide_legend(override.aes = list(linetype = c("solid", "solid"))))

###########################