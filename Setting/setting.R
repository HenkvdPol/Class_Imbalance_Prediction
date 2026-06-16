################################################################################
# Simulation settings
################################################################################

# Levels of sample sizes
n <- c(3000, 10^4)

# levels of prevalences
prevalence <- c(0.025, 0.05, 0.10, 0.20)

# Variance of the covariates
sigma <- 1

# Levels of interaction between the covariates
delta <- 0 #c(0, 0.25, 0.5)

# Levels of dimension
dimension <- c(1,2,4,8,16)

# Number of Monte carlo repeats.
R <- 250

# Use parallel computing or not.
use_cluster <- TRUE

# The measures that are computed for prediction.
Measures <- c("Accuracy", "Precision", "Recall","Specificity", 
              "FNR","FPR","F1",  "AUC",  "APR")

# The models that are used for prediction.
Models <- c('Original', 'ROS', 'RUS', 'SMOTE', 'Prevalence-threshold',
            'AUC-threshold', 'F1-threshold', 'APR-threshold')

# Starting seed for reproducibility
startSeed0 <- 111111

# The parameter Beta. We set 8 coefficients to 0 to introduce some noise.
coef <- c(0, rnorm(15, 0, 0.5))
coef[c(3,7,11,15)] <- 0

# The scenarios. Here; Beta_1 = coef[1] is increased by 0.2 every scenario.
Scenario <- lapply(2 * (0:10) /10, function(x){coef + c(coef[1] + x, 
                                                        rep(0,length(coef) - 1))})
