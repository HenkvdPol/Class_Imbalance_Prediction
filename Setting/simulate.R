################################################################################
# Simulation settings
################################################################################

set.seed(1)
source('setting.R')

if(length(commandArgs(trailingOnly = TRUE)) == 0){
  # Test setting when in R and don't want to use the cluster.
  args <- c(3000, 2, 0, 12, TRUE, 1)
  use_cluster <- as.logical(args[5])
  R <- as.integer(args[6])
}else{
  args <- commandArgs(trailingOnly = TRUE)
}

cat(args)

n <- as.integer(args[1]) # sample size
dimension <- as.integer(args[2]) # dimension
delta <- as.double(args[3])
extra_seed <- as.integer(args[4])
use_cluster <- as.logical(args[5])
R <- as.integer(args[6])

cat('\n', n, dimension, delta, extra_seed)

sim <- expand.grid('prevalence' = prevalence,
                   'delta' = delta,
                   'sample_size' = n,
                   'Scenario' = Scenario,
                   'dimension' = dimension)

sim$id <- 1:nrow(sim)

################################################################################
# Libraries and additional functions to create simulation.
################################################################################

# Necessary libraries
library(ROCR)
library(caret)
library(PRROC)
library(pROC)
library(dbplyr)
library(foreach)
library(doParallel)
library(smotefamily)
library(CalibrationCurves)

# Additional functions in use.
source('additional_functions.R')


################################################################################
# Simulation function 
################################################################################

genSim <- function(R, n, prev, Beta, startSeed, p, delta, sigma, id){
  # R = number of Monte Carlo repeats.
  # n = sample size. Here either 3000 or 10^4
  # prev = prevalence level in [0,1]
  # Beta = coefficients
  # startSeed = seed for reproducibility
  # p = covariate dimension
  # delta = interaction parameter of the covariates, in [0,1]
  # sigma = variance of the covariates
  # id = the simulation number
  cat('\n seed number;', startSeed)
  
  
  # Beta depends on the dimension
  Beta <- Beta[1:p]
  cat('\n Beta = ', Beta)

  ## estimated intercept for the desired prevalence
  b0 <- targetProp.fun(prev = prev, Beta = Beta, setSeed = startSeed, 
                        p = p, delta = delta, sigma = sigma)
  cat(paste('\n intercept =', b0))
  
  ## generate the covariate matrix to be keep fixed 
  X <- generate_X(n = n, p = p, delta = delta, sigma = 1, setSeed = startSeed)
  
  # Initialize the individual outcome arrays. Note the order.
  Result <- expand.grid('Measures' = Measures,
                        'Model' = Models,
                        'R' = 1:R)
  Result$Value <- NA
  Result$Calibration <- NA
  
  # We also add the parameters here to store the information.
  Result$Prevalence <- prev
  Result$Dimension <- p
  Result$delta <- delta
  Result$sigma <- sigma
  Result$n <- n
  Result$Scenario <- which(vapply(Scenario, function(x) all(x[1:p] == Beta), logical(1)))
  Result$Simulation <- id
  
  
  # Also collect the true mean of the sample and the selected threshold.
  Result$`True prevalence` <- NA
  Result$Threshold <- 1/2 # it is either 1/2, or something else.
  
  #######################################################################
  ################### Start Simulation  #################################
  #######################################################################
  
  for (i in 1:R) {
    cat('\n id =', id, '\n', 'i =', i)
    
    # Set seed
    seed_i <- 111111 + startSeed + i
    
    # generate outcome and error-handling
    
    y <- generate_outcome(Beta0 = b0, 
                          X = X, 
                          BetasX = Beta, 
                          setSeed = seed_i)
    
    # Form the dataset (X,Y)
    dataComplete <- as.data.frame(cbind(X, y))
    names(dataComplete) <- c('Intercept', paste('x', 1:(dim(X)[2] - 1), sep = ''), 'y')
    form <- as.formula(paste0('y ~ ', paste('x', 1:(dim(X)[2] - 1), sep = '', collapse = ' + ')))
    
    # Split the dataset depending on the outcome.
    indDiv <- createDataPartition(dataComplete$y, p =  2/3, list = FALSE)
    datatrain <- dataComplete[indDiv, ]
    datatest  <- dataComplete[-indDiv, ]
    
    #print('hello1')
    
    
    # Error-handling
    l <- 0
    while(sum(y) < 10 || sum(datatrain$y) < 5 || sum(datatest$y) < 2){
      seed_i <- seed_i + 1
      
      # If we don't find a new one, we stop. Then it is not possible.
      # We also make here note of this. But this does not happen.
      if(l > 50){
        
        cat('\n Not enough events!')
        
        return(Result)
      }
      
      cat(paste('\n error-handling', l))
      
      l <- l + 1
      
      y <- generate_outcome(Beta0 = b0, 
                            X = X, 
                            BetasX = Beta, 
                            setSeed = seed_i)
      
      # Form the dataset (X,Y)
      dataComplete <- as.data.frame(cbind(X, y))
      names(dataComplete) <- c('Intercept', paste('x', 1:(dim(X)[2] - 1), sep = ''), 'y')
      
      #print(head(dataComplete))
      
      # Split the dataset depending on the outcome.
      indDiv <- createDataPartition(dataComplete$y, p =  0.68, list = FALSE )
      datatrain <- dataComplete[indDiv, ]
      datatest  <- dataComplete[-indDiv, ]
    }
    rm(l)
    
    # Collect the observed prevalence.
    k <- length(Models) * length(Measures)
    
    # Here we use the order of the dataframe.
    Result[(i-1)*k + (1:k), 'True prevalence'] <- mean(y)
    cat(paste('\n observed mean(y) =', mean(y)))
    
    #print(paste('mean outcome train:', mean(datatrain$y),', mean outcome test:', mean(datatest$y)))
    #print(paste('outcome y train:', sum(datatrain$y),', outcome y test:', sum(datatest$y)))
    
    #print('hello2')
    
    # Re-sampling
    
    # ROS dataset
    trainup <- upSample(x = datatrain[, -dim(datatrain)[2]], y = as.factor(datatrain$y))
    trainup$y <-  as.numeric(as.character(trainup$Class))
    trainup <- trainup[,names(trainup) !="Class" ]
    
    #print('hello3')
    
    # RUS dataset
    traindown <- downSample(x = datatrain[, -dim(datatrain)[2]], y = as.factor(datatrain$y))
    traindown$y <-  as.numeric(as.character(traindown$Class))
    traindown <- traindown[,names(traindown) !="Class" ]
    
    #print('hello4')
    
    # SMOTE dataset
    event_index <- which(datatrain$y == 1)
    y_index <- which(names(datatrain) == 'y')
    trainsmote <- datatrain
    
    for(j in 1:(sum(datatrain$y == 0) - sum(datatrain$y == 1))){
      new_event <- colSums(datatrain[sample(event_index, 5), 
                                     -y_index]) / 5
      trainsmote <- rbind(trainsmote, c(new_event, 1))
    }
    
    #print('hello5')
    
    # Prediction
    pred_train <- fit_lr(datatrain, datatest)
    pred_trainup <- fit_lr(trainup, datatest)
    pred_traindown <- fit_lr(traindown, datatest)
    pred_trainsmote <- fit_lr(trainsmote, datatest)
    
    #print('hello6')
    
    
    # Compute performance measures. The models are in order;
    # (1) Original, (2) ROS, (3) RUS, (4) SMOTE, (5) Prevalence-threshold
    # (6) AUC-threshold, (7) F1-threshold, (8) APR-threshold
    m <- length(Measures)
    Result$Value[(i-1)*k + ((m*1 - m) +1):(m*1)] <- performanceMeasures(pred_train)
    Result$Calibration[(i-1)*k + ((m*1 - m) +1):(m*1)] <- compute_cal(pred_train)
    
    Result$Value[(i-1)*k + ((m*2 - m) +1):(m*2)] <- performanceMeasures(pred_trainup)
    Result$Calibration[(i-1)*k + ((m*2 - m) +1):(m*2)] <- compute_cal(pred_trainup)
    
    Result$Value[(i-1)*k + ((m*3 - m) +1):(m*3)] <- performanceMeasures(pred_traindown)
    Result$Calibration[(i-1)*k + ((m*3 - m) +1):(m*3)] <- compute_cal(pred_traindown)
    
    Result$Value[(i-1)*k + ((m*4 - m) +1):(m*4)] <- performanceMeasures(pred_trainsmote)
    Result$Calibration[(i-1)*k + ((m*4 - m) +1):(m*4)] <- compute_cal(pred_trainsmote)
    
    p#print('hello7')
    
    Result$Value[(i-1)*k + ((m*5 - m) +1):(m*5)] <- performanceMeasures(pred_train,
                                                                        mean(datatrain$y))
    Result$Threshold[(i-1)*k + ((m*5 - m) +1):(m*5)] <- mean(datatrain$y)
    
    #print('hello8')
    
    auc_threshold <- threshold_auc(pred_train$prediction, datatest$y)
    Result$Value[(i-1)*k + ((m*6 - m) +1):(m*6)] <- performanceMeasures(pred_train,
                                                                        auc_threshold)
    Result$Threshold[(i-1)*k + ((m*6 - m) +1):(m*6)] <- auc_threshold
    
    #print('hello9')
    
    f1_threshold <- threshold_f1(pred_train$prediction, datatest$y)
    Result$Value[(i-1)*k + ((m*7 - m) +1):(m*7)] <- performanceMeasures(pred_train,
                                                                        f1_threshold)
    Result$Threshold[(i-1)*k + ((m*7 - m) +1):(m*7)] <- f1_threshold
    
    #print('hello10')
    
    apr_threshold <- threshold_apr(pred_train$prediction, datatest$y)
    Result$Value[(i-1)*k + ((m*8 - m) +1):(m*8)] <- performanceMeasures(pred_train,
                                                                        apr_threshold)
    Result$Threshold[(i-1)*k + ((m*8 - m) +1):(m*8)] <- apr_threshold
    
    #print('hello11')
  }
  
  # We save the results for every iteration.
  save(Result, file = paste0('Result/sim-',n ,'-',p,'-',extra_seed,'-',R,'-',id, '.Rdata'))
  
  # Return the result.
  return(id)
}

if(use_cluster){
message("SLURM_CPUS_PER_TASK: ", Sys.getenv("SLURM_CPUS_PER_TASK"))
message("SLURM_JOB_ID: ", Sys.getenv("SLURM_JOB_ID"))
message("detectCores: ", parallel::detectCores())

  # Want to know the time
  start <- proc.time()
  # Make clusters to parallel
  n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
  registerDoParallel(cores = n_cores)  # no makeCluster at all
  
  message("Registered parallel backend: ", getDoParName())
  message("Number of workers: ", getDoParWorkers())
  
  tmp <- list()
  
  foreach(i=1:nrow(sim), 
                 .packages = c('ROCR', 'caret', 'PRROC', 'pROC', 'MASS', 
                               'dplyr', 'CalibrationCurves')) %dopar% {genSim(R, 
                                                         sim$sample_size[i], 
                                                         sim$prevalence[i], 
                                                         sim$Scenario[[i]], 
                                                         i + extra_seed, # the start seed is the id. 
                                                         sim$dimension[i], 
                                                         sim$delta[i], 
                                                         sigma,
                                                         i)}
  
  # Stop clusters.
  stopImplicitCluster() 
  print(proc.time() - start)
}

if(!use_cluster){
  start <- proc.time()
  
  tmp <- list()
  
  for(i in 1:nrow(sim)){
    genSim(R,
           sim$sample_size[i],
           sim$prevalence[i],  
           sim$Scenario[[i]],
           i + extra_seed, # The start seed is the id 
           sim$dimension[i], 
           sim$delta[i],
           sigma,
           i)
  }
}
