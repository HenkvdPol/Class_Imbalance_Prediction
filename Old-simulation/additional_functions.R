# File with additional functions


# Necessary libraries
library(ROCR)
library(caret)
library(brglm2)
library(RSBID)
library(PRROC)
library(pROC)
library(dbplyr)
library(foreach)
library(doParallel)
library(smotefamily)


# Function to generate the covariate matrix
generate_X <- function(n, setSeed = 123, old = FALSE){
  set.seed(setSeed)
  
  x1 <- rnorm(n, 0, 1)
  x2 <- rnorm(n, 0, 1)
  x3 <- rpois(n, 4)
  x4 <- runif(n, 0, 1)
  x5 <- rnorm(n, 0, 1)
  x6 <- rexp(n, 3)
  x7 <- rnorm(n, 0, 1)
  x8 <- rexp(n,1/4) # x8p <- 1/x8
  x9 <- rnorm(n, 0, 1)
  
  
  # May 2025; add the standardized results.
    x1 <- scale(x1)
    x2 <- scale(x2)
    x3 <- scale(x3)
    x4 <- scale(x4)
    x5 <- scale(x5)
    x6 <- scale(x6)
    x7 <- scale(x7)
    x8 <- scale(1/x8)
    x9 <- scale(x9)
  
  X <- model.matrix(~ x1+x2+x3+x4+x5+x6+x7+x8+x9)
  
  if(old){
    x1 <- rnorm(n, 0, 1)
    x2 <- rnorm(n, 0, 1)
    x3 <- as.factor(rbinom(n, 1, 0.5))
    x4 <- runif(n, 0, 1)
    x5 <- rnorm(n, 0, 1)
    x6 <- as.factor(rbinom(n, 1, 0.5))
    x7 <- rnorm(n, 0, 1)
    x8 <- rexp(n,1/4)
    x9 <- as.factor(sample(1:3, size = n, replace = T))
    X <- model.matrix(~ x1+x2+x3+x4+x5+x6+x7+I(1/x8)+x9+I(1/x8):x9)
  }
  
  return(X)
}



# Compute the intercept.
targetProp.fun <- function(prop, param,setSeed, old = FALSE){
  N.pop <- 1000000
  # Size of super-population.
  
  # Target prevalence of the outcome in the super-population.
  target.prop.outcome <- prop
  
  # How close is the empirical outcome prevalence required to be to the
  # target outcome prevalence.
  tolerance <- 0.0001
  
  ################################################################################
  # Generate 9 baseline covariates for each subject in the super-population.
  # The 7 are from independent standard normal distibutions.
  ################################################################################
  
  X <- generate_X(N.pop, setSeed = setSeed, old = old)
  
  ################################################################################
  # Generate a binary treatment variable with the given prevalence of treatment.
  ################################################################################

  B.outcome <- param
  
  ################################################################################
  # Define a function for generating binary outcomes with a given value of
  # the intercept in the outcomes model.
  ################################################################################
  
  outcome.function <- function(b0.outcome){
    
    beta.outcome.modified <- c(b0.outcome,B.outcome)
    # Set the intercept of the outcome model to the given value.
    
    XB <- X %*% beta.outcome.modified
    # Linear predictor.
    
    p.outcome <- exp(XB)/(1 + exp(XB))
    # Probability of the outcome.
    
    Y <- rbinom(N.pop,1,p.outcome)
    
    return(Y)
    
    remove(beta.outcome.modified,XB,p.outcome,Y)
    
  }
  
  ################################################################################
  # Use a bisection approach to determine the intercept that results in the
  # desired prevalence of the outcome.
  ################################################################################
  
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
    
    outcome <- outcome.function(b0.outcome=int.mid)
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


#########################################################################
###### find threshold that maximize both sensitivity and specificity ####
#########################################################################
threshold <- function(predict, response) {
  r <- pROC::roc(response, predict)
  d <- sqrt((r$sensitivities-1)^2+(1-r$specificities)^2)
  r$thresholds[which.min(d)]
}

#########################################################################
##### find threshold that maximize F1 measure                     #######
#########################################################################

threshold_f1 <- function(predict, response) {
  thvalues <- seq(0.001,1-0.001,by=0.001)
  f1values <- sapply(thvalues, function(t) {
    out <-confusionMatrix(data=as.factor(as.numeric(predict > t )), 
                          reference = factor(response),
                          positive = "1")
    out$byClass[c("Precision","Recall","F1")]
  })
  f1values <- as.data.frame(t(f1values))
  res <- thvalues[which.max(f1values$F1)]
  if(length(res) == 0){res <- 0}
  
  return(res)
} 

#######################################################################
####### find threshold for which both precision and recall are high ###
#######################################################################

threshold_pr <- function(predict, response) {
  thvalues <- seq(0.001,1-0.001,by=0.001)
  f1values <- sapply(thvalues, function(t) {
    out <-confusionMatrix(data= factor(as.numeric(predict > t ),
                                       levels = c(0,1)), 
                          reference = factor(response, levels = c(0,1)),
                          positive = "1")
    out$byClass[c("Precision","Recall","F1")]
  })
  f1values <- as.data.frame(t(f1values))
  res <- apply(f1values[,-3],1,function(t) sqrt(sum((t-c(1,1))^2)))
  th =thvalues[which.min(res )]
} 

#######################################################################
####### find threshold for one that maximizes the APR curve ###########
#######################################################################

threshold_APR_curve <- function(predict, response) {
  pr <- pr.curve(predict[response == 1], 
                 predict[response == 0], 
                 curve = T)
  
  best_pr <- apply(pr$curve[,-3],1,function(t) sqrt(sum((t-c(1,1))^2)))
  threshold <- pr$curve[,3][which.min(best_pr )]
  
  return(threshold)
} 


## generate outcome
genOutcomey <- function(beta0, X, betasX, setSeed) {
  set.seed(setSeed)
  n <- nrow(X)
  beta.modified <- c(beta0, betasX)
  ### linear predictor.
  XA <- X %*% beta.modified
  ### probability of the treatment.
  p.treat <- exp(XA)/(1 + exp(XA))
  ### generation the treatment groups
  out <- rbinom(n, 1, p.treat)
  return(out)
}

#### fit logistic regression model #####
modreg <- function(trainingData, testingData) {
  fit <- glm(formula = y~., family = 'binomial', data = trainingData)
  #print(summary(fit))
  pred <- predict( fit, newdata = testingData, type = "response")
  testingData$prediction <- pred 
  return(testingData)
} 

#### performance measures function #####
performanceMeasures <- function(datapred, threshold = 0.5) {
  roc <- unlist(slot(performance(prediction(datapred$prediction, datapred$y), "auc"), "y.values"))
  
  out <-confusionMatrix(data= factor(as.numeric(datapred$prediction > threshold),
                                       levels = c(0,1)), 
                        reference = factor(datapred$y, levels = c(0,1)),
                        positive = "1",
                        mode = "everything")
  
  apr <- pr.curve(scores.class0 = unlist(datapred$prediction), 
                  weights.class0 = datapred$y, 
                  curve = FALSE)$auc.integral
  
  
  out <- c(out$overall[1], out$byClass["Precision"], out$byClass["Recall"],
           out$byClass["Specificity"], 1 - out$byClass["Recall"],
           1 - out$byClass["Specificity"], out$byClass["F1"], roc, apr)
  names(out)[c(5,6,7,8,9)] <-c("FNR","FPR","F1","AUC", "APR")
  
  return(out)
}
