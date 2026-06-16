# File with additional functions


# Necessary libraries
library(ROCR)
library(caret)
library(PRROC)
library(pROC)
library(dbplyr)
library(foreach)
library(doParallel)
library(MASS)
library(Matrix)
library(CalibrationCurves)

################################################################################
# Function to generate the covariate matrix
################################################################################
generate_X <- function(n, p, delta, sigma, setSeed){
  set.seed(setSeed)
  
  cat(paste('\n p in generate x;', p))
  
  if(p == 1){
    X <- rnorm(n=n, mean = 0, sd = sigma)
    X <- cbind(1, X)
    return(X)
  }
  
  mu <- rep(0, p) # mean vector
  delta <- abs(delta) # delta needs to be positive. For safety
  block <- matrix(c(sigma, delta * sigma,
                          delta * sigma, sigma),
                        nrow = 2, byrow = TRUE)
  
  #cat('\n', block)
  
  Sigma <- kronecker(diag(as.integer(p/2)), block)
  
  #cat('\n', Sigma)
  
  X <- mvrnorm(n =n, mu = mu, Sigma = Sigma)
  X <- cbind(1, X)

  return(X)
}


################################################################################
# Compute the intercept.
################################################################################
targetProp.fun <- function(prev, Beta, setSeed, p, delta, sigma){
  N.pop <- 1000000
  # Size of super-population.
  
  # Target prevalence of the outcome in the super-population.
  target.prop.outcome <- prev
  
  # How close is the empirical outcome prevalence required to be to the
  # target outcome prevalence.
  tolerance <- 0.0001
  
  ##############################################################################
  # Generate 9 baseline covariates for each subject in the super-population.
  # The 7 are from independent standard normal distibutions.
  ##############################################################################
  
  X <- generate_X(n = N.pop, setSeed = setSeed, 
                  p = p, delta = delta, sigma = sigma)
  
  ##############################################################################
  # Define a function for generating binary outcomes with a given value of
  # the intercept in the outcomes model.
  ##############################################################################
  
  outcome.function <- function(b0.outcome, Beta = Beta){
    
    Beta.outcome.modified <- c(b0.outcome, Beta)
    # Set the intercept of the outcome model to the given value.
    
    XB <- X %*% Beta.outcome.modified
    # Linear predictor.
    
    p.outcome <- exp(XB)/(1 + exp(XB))
    # Probability of the outcome.
    
    Y <- rbinom(N.pop,1,p.outcome)
    
    return(Y)
    
    remove(Beta.outcome.modified,XB,p.outcome,Y)
    
  }
  
  ##############################################################################
  # Use a bisection approach to determine the intercept that results in the
  # desired prevalence of the outcome.
  ##############################################################################
  
  # Define endpoints of interval to bisect.
  int.low <- -30
  int.high <- 30
  
  iter <- 1
  
  outcome.prev <- 1
  # Initial value of the empirical value of outcome prevalence.
  # This can be any value as long as it differs from the target value.
  
  while(abs(outcome.prev - target.prop.outcome) > tolerance){
    set.seed(iter)
    int.mid <- (int.low + int.high)/2
    
    outcome <- outcome.function(b0.outcome=int.mid, Beta = Beta)
    outcome.prev <- mean(outcome)
    
    if (outcome.prev < target.prop.outcome) int.low <- int.mid else
      int.high <- int.mid
    
    iter <- iter + 1
    
  }
  
  return(int.mid)
  # Intercept for the regression model to produce the desired outcome prevalence.
  
  remove(int.low,int.mid,int.high,iter,outcome.prev,outcome.function,
         B.outcome,X)
}


################################################################################
# Find threshold that maximize both sensitivity and specificity (AUC) 
################################################################################
threshold_auc <- function(predict, response) {
  r <- pROC::roc(response, predict)
  d <- sqrt((r$sensitivities-1)^2+(1-r$specificities)^2)
  r$thresholds[which.min(d)]
}

################################################################################
# Find threshold that maximize F1 measure                
################################################################################

threshold_f1 <- function(predict, response) {
  thvalues <- seq(0.001,1-0.001,by=0.001)
  f1values <- sapply(thvalues, function(t) {
    out <-suppressWarnings(
      confusionMatrix(data=as.factor(as.numeric(predict > t )), 
                      reference = factor(response),
                      positive = "1"))
    out$byClass[c("Precision","Recall","F1")]
  })
  f1values <- as.data.frame(t(f1values))
  res <- thvalues[which.max(f1values$F1)]
  if(length(res) == 0){res <- 0}
  
  return(res)
} 

################################################################################
# Find threshold for which both precision and recall are high
################################################################################

threshold_pr <- function(predict, response) {
  thvalues <- seq(0.001,1-0.001,by=0.001)
  f1values <- sapply(thvalues, function(t) {
    out <-suppressWarnings(
      confusionMatrix(data= factor(as.numeric(predict > t ),
                                   levels = c(0,1)), 
                      reference = factor(response, levels = c(0,1)),
                      positive = "1"))
    out$byClass[c("Precision","Recall","F1")]
  })
  f1values <- as.data.frame(t(f1values))
  res <- apply(f1values[,-3],1,function(t) sqrt(sum((t-c(1,1))^2)))
  th =thvalues[which.min(res )]
} 

################################################################################
# Find threshold for one that maximizes the APR curve 
################################################################################

threshold_apr <- function(predict, response) {
  pr <- pr.curve(predict[response == 1], 
                 predict[response == 0], 
                 curve = T)
  
  best_pr <- apply(pr$curve[,-3],1,function(t) sqrt(sum((t-c(1,1))^2)))
  threshold <- pr$curve[,3][which.min(best_pr )]
  
  return(threshold)
} 

################################################################################
# generate outcome using the logistic regression model
################################################################################
generate_outcome <- function(Beta0, X, BetasX, setSeed){
  set.seed(setSeed)
  n <- nrow(X)
  Beta.modified <- c(Beta0, BetasX)
  ### linear predictor.
  XB <- X %*% Beta.modified
  ### probability of the treatment.
  p.treat <- exp(XB)/(1 + exp(XB))
  ### generation the treatment groups
  out <- rbinom(n, 1, p.treat)
  return(out)
}

################################################################################
# fit data to the logistic regression model
################################################################################
fit_lr <- function(trainingData, testingData) {
  fit <- suppressWarnings(glm(formula = y~., family = 'binomial', data = trainingData))
  #print(summary(fit))
  
  pred <- predict( fit, newdata = testingData, type = "response")
  testingData$prediction <- pred
  
  # Also check if it converged.
  testingData$converged <- fit$converged
    
  return(testingData)
} 

################################################################################
# performance measures function
################################################################################
performanceMeasures <- function(datapred, threshold = 0.5) {
  # AUC
  roc <- unlist(slot(performance(prediction(datapred$prediction, datapred$y),
                                 "auc"), "y.values"))
  
  # FNR, FPR, F1
  out <-suppressWarnings(
    confusionMatrix(data= factor(as.numeric(datapred$prediction > threshold),
                                 levels = c(0,1)),
                    reference = factor(datapred$y, levels = c(0,1)),
                    positive = "1",
                    mode = "everything"))
  
  # APR
  apr <- pr.curve(scores.class0 = unlist(datapred$prediction), 
                  weights.class0 = datapred$y, 
                  curve = FALSE)$auc.integral
  
  # set everything in one array.
  out <- c(out$overall[1], out$byClass["Precision"], out$byClass["Recall"],
           out$byClass["Specificity"], 1 - out$byClass["Recall"],
           1 - out$byClass["Specificity"], out$byClass["F1"], roc, apr)
  names(out)[c(5,6,7,8,9)] <-c("FNR","FPR","F1","AUC", "APR")
  
  return(out)
}

################################################################################
# compute the calibration              
################################################################################

# Since the fitted probabilities may be general prone to error, we adjust here
# the val.prob.ci.2 function here to have degree = 1. However, this is also 
# needed for one of the wrapper functions BT.samples.
#vpc2 <- val.prob.ci.2
#
#BT.samples <- function(y, p, to.pred){
#  Df = cbind.data.frame(y, p)
#  repeat {
#    BT.sample = Df[sample(1:nrow(Df), replace = T), ]
#    loess.BT = loess(y ~ p, BT.sample, degree = 1) # here we add degree = 1.
#    pred.loess = predict(loess.BT, to.pred, type = "fitted")
#    if (!any(is.na(pred.loess))) 
#      break
#  }
#  return(pred.loess)
#}



compute_cal <- function(datapred){
  
  if(!all(datapred$converged)){
    # If the algorithm did not converge, we can not compute the calibration
    return(NA)
  }
  
  #print('hello-a')
  
  # Check if we receive non-NAN fitted values
  fit <- suppressWarnings(loess(y ~ prediction, data = datapred))
  if(anyNA(fit$fitted)){
    return(NA)
  }
  
  #print('hello-b')
  
  # This can also happen in the bootstrap fitting.
  to.pred <- seq(min(datapred$prediction), 
                 max(datapred$prediction), 
                 length = 200)
  # How BT.loess is working
  Df = cbind.data.frame(y = datapred$y, p = datapred$prediction)
  BT.sample = Df[sample(1:nrow(Df), replace = T), ]
  loess.BT = suppressWarnings(loess(y ~ p, BT.sample))
  if(anyNA(loess.BT$fitted)){
    return(NA)
  }
  
  #print('hello-c')
  
  # Or just prediction goes wrong in general
  SmFit = suppressWarnings(loess(y ~ prediction, data = datapred))
  if(anyNA(SmFit$fitted)){
    return(NA)
  }
  #cl.loess = predict(SmFit, type = "fitted", se = TRUE)
  
  rm(fit, to.pred, Df, BT.sample, loess.BT, SmFit)
  
  #print('hello-d')
  
  
  # Otherwise, we continue as usual. We cap the probabilities.
  eps <- 1e-8
  datapred$prediction <- pmin(pmax(datapred$prediction, eps), 1 - eps)
  
  # Since val.prob.ci.2() will only return calibration curves if we set 
  # pl = TRUE, we need to create the png and then immedeatly delete it, so
  # we do not clutter the memory
  pdf(NULL)
  cal <- val.prob.ci.2(datapred$prediction, datapred$y, pl = TRUE)
  dev.off()
  
  # We round the numbers down to 4 digits, as we are only interested in the
  # plot
  horizontal <- round(cal$CalibrationCurves$FlexibleCalibration$x,8)
  horizontal <- paste(horizontal, collapse = ', ')
  
  vertical <- round(cal$CalibrationCurves$FlexibleCalibration$y, 8)
  vertical <- paste(vertical, collapse = ', ')
  
  # We only collect the point estimate here.
  intercept <- unname(cal$Calibration$Intercept[1])
  slope <- unname(cal$Calibration$Slope[1])
  
  out <- c('x', horizontal, 'y', vertical, 'intercept', intercept,
           'slope', slope, 'check')
  
  rm(intercept, slope, horizontal, vertical, cal)
  return(out)
}

################################################################################
# test gen simulation
################################################################################

genSim_test <- function(R, n, prev, Beta, startSeed, p, delta, sigma, id){
  # R = number of Monte Carlo repeats.
  # n = sample size. Here either 10^3 or 10^4
  # prev = prevalence level in [0,1]
  # Beta = coefficients
  # startSeed = seed for reproducibility
  # p = covariate dimension
  # delta = interaction parameter of the covariates, in [0,1]
  # sigma = variance of the covariates
  # id = the simulation number
  
  
  # Beta depends on the dimension
  Beta <- Beta[1:p]
  
  ## estimated intercept for the desired prevalence
  b0 <- targetProp.fun(prev = prev, Beta = Beta, setSeed = startSeed, 
                       p = p, delta = delta, sigma = sigma)
  cat(print(paste('intercept', b0)))
  
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
    #cat(paste('id', id))
    #cat(paste('i', i))
    
    # Set seed
    seed_i <- 111111 + i
    
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
    
    print('hello1')
    
    
    # Error-handling
    l <- 0
    while(sum(y) < 10 || sum(datatrain$y) < 5 || sum(datatest$y) < 2){
      seed_i <- seed_i + 1
      
      # If we don't find a new one, we stop. Then it is not possible.
      if(l > 1000){return(Result)}
      
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
    
    # Collect the observed prevalence.
    k <- length(Models) * length(Measures)
    
    # Here we use the order of the dataframe.
    Result[(i-1)*k + (1:k), 'True prevalence'] <- mean(y)
    print(mean(y))
    
    print(paste('mean outcome train:', mean(datatrain$y),
                ', mean outcome test:', mean(datatest$y)))
    print(paste('outcome y train:', sum(datatrain$y),
                ', outcome y test:', sum(datatest$y)))
    
    print('hello2')
    
    # Re-sampling
    
    # ROS dataset
    trainup <- upSample(x = datatrain[, -dim(datatrain)[2]], y = as.factor(datatrain$y))
    trainup$y <-  as.numeric(as.character(trainup$Class))
    trainup <- trainup[,names(trainup) !="Class" ]
    
    print('hello3')
    
    # RUS dataset
    traindown <- downSample(x = datatrain[, -dim(datatrain)[2]], y = as.factor(datatrain$y))
    traindown$y <-  as.numeric(as.character(traindown$Class))
    traindown <- traindown[,names(traindown) !="Class" ]
    
    print('hello4')
    
    # SMOTE dataset
    event_index <- which(datatrain$y == 1)
    y_index <- which(names(datatrain) == 'y')
    trainsmote <- datatrain
    
    for(j in 1:(sum(datatrain$y == 0) - sum(datatrain$y == 1))){
      new_event <- colSums(datatrain[sample(event_index, 5), 
                                     -y_index]) / 5
      trainsmote <- rbind(trainsmote, c(new_event, 1))
    }
    
    print('hello5')
    
    # Prediction
    pred_train <- fit_lr(datatrain, datatest)
    pred_trainup <- fit_lr(trainup, datatest)
    pred_traindown <- fit_lr(traindown, datatest)
    pred_trainsmote <- fit_lr(trainsmote, datatest)
    
    print('hello6')
    
    
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
    
    print('hello7')
    
    Result$Value[(i-1)*k + ((m*5 - m) +1):(m*5)] <- performanceMeasures(pred_train,
                                                                        mean(datatrain$y))
    Result$Threshold[(i-1)*k + ((m*5 - m) +1):(m*5)] <- mean(datatrain$y)
    
    print('hello8')
    
    auc_threshold <- threshold_auc(pred_train$prediction, datatest$y)
    Result$Value[(i-1)*k + ((m*6 - m) +1):(m*6)] <- performanceMeasures(pred_train,
                                                                        auc_threshold)
    Result$Threshold[(i-1)*k + ((m*6 - m) +1):(m*6)] <- auc_threshold
    
    print('hello9')
    
    f1_threshold <- threshold_f1(pred_train$prediction, datatest$y)
    Result$Value[(i-1)*k + ((m*7 - m) +1):(m*7)] <- performanceMeasures(pred_train,
                                                                        f1_threshold)
    Result$Threshold[(i-1)*k + ((m*7 - m) +1):(m*7)] <- f1_threshold
    
    print('hello10')
    
    apr_threshold <- threshold_apr(pred_train$prediction, datatest$y)
    Result$Value[(i-1)*k + ((m*8 - m) +1):(m*8)] <- performanceMeasures(pred_train,
                                                                        apr_threshold)
    Result$Threshold[(i-1)*k + ((m*8 - m) +1):(m*8)] <- auc_threshold
    
    print('hello11')
  }
  
  # Return the result.
  return(Result)
}