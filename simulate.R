# Sim settings.
n <- 10000
prevalence <- c(0.025, 0.05, 0.10, 0.20)
R <- 1000
use_cluster <- TRUE


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

# Additional functions in use.
source('additional_functions.R')

#######################################################################
################### Simulation function  ##############################
#######################################################################

genSim <- function(R, n, prev, betas_x, startSeed){

  ## estimated intercept for the desired prevalence 
  seed <- startSeed
  b0 <- targetProp.fun(prop = prev, param = betas_x, setSeed = seed)
  
  ## generate the covariate matrix to be keep fixed 
  X <- generate_X(n, setSeed = seed)
  
  # Initialize the individual outcome arrays.
  olr <- ros <- rus <- smt <- NULL
  olr_br <- ros_br <- rus_br <- smt_br <- NULL
  olr_corrected <-olr_bth_sp <- olr_bth_f1 <-NULL
  
  olr_bth_APR_curve <- NULL
  
  meanY <- NULL
  
  #######################################################################
  ################### Start Simulation  #################################
  #######################################################################
  
  for (i in 1:R) {
    seed_i <- 111111 + i
    
    # generate outcome and error-handling
    
    y <- genOutcomey(beta0 = b0, 
                     X = X, 
                     betasX = betas_x, 
                     setSeed = seed_i)
    
    # Form the dataset (X,Y)
    dataComplete <- as.data.frame(cbind(X, y))
    names(dataComplete) <- c('Intercept', paste('x', 1:(dim(X)[2] - 1), sep = ''), 'y')
    
    # Split the dataset depending on the outcome.
    indDiv <- createDataPartition(dataComplete$y, p =  0.68, list = FALSE )
    datatrain <- dataComplete[indDiv, ]
    datatest  <- dataComplete[-indDiv, ]
    
    
    # Error-handling
    while(sum(y) < 10 || sum(datatrain$y) < 5 || sum(datatest$y) < 2){
      seed_i <- seed_i + 1
      print('hello')
      
      y <- genOutcomey(beta0 = b0, 
                       X = X, 
                       betasX = betas_x, 
                       setSeed = seed_i)
      
      # Form the dataset (X,Y)
      dataComplete <- as.data.frame(cbind(X, y))
      names(dataComplete) <- c('Intercept', paste('x', 1:(dim(X)[2] - 1), sep = ''), 'y')
      
      # Split the dataset depending on the outcome.
      indDiv <- createDataPartition(dataComplete$y, p =  0.68, list = FALSE )
      datatrain <- dataComplete[indDiv, ]
      datatest  <- dataComplete[-indDiv, ]
    }
    
    meanY <- c(meanY, mean(y))
    print(mean(y))
    
    print(paste('mean outcome train:', mean(datatrain$y),
                ', mean outcome test:', mean(datatest$y)))
    print(paste('outcome y train:', sum(datatrain$y),
           ', outcome y test:', sum(datatest$y)))
    
    # Re-sampling
    
    # ROS dataset
    trainup <- upSample(x = datatrain[, -dim(datatrain)[2]], y = as.factor(datatrain$y))
    trainup$y <-  as.numeric(as.character(trainup$Class))
    trainup <- trainup[,names(trainup) !="Class" ]
    
    # RUS dataset
    traindown <- downSample(x = datatrain[, -dim(datatrain)[2]], y = as.factor(datatrain$y))
    traindown$y <-  as.numeric(as.character(traindown$Class))
    traindown <- traindown[,names(traindown) !="Class" ]
    
    # SMOTE dataset
    k <- ifelse(sum(datatrain$y) > 5, 5, sum(datatrain$y) - 1)
    
    trainsmote <- RSBID::SMOTE(datatrain, 'y', perc_maj=100, k=k)
    
    # Prediction

    pred_train <- modreg(datatrain, datatest)
    pred_trainup <- modreg(trainup, datatest)
    pred_traindown <- modreg(traindown, datatest)
    pred_trainsmote <- modreg(trainsmote, datatest)
    
    
    # Compute performance measures
    
    olr <- rbind(olr, performanceMeasures(pred_train))
    ros <- rbind(ros, performanceMeasures(pred_trainup))
    rus <- rbind(rus, performanceMeasures(pred_traindown))
    smt <-  rbind(smt, performanceMeasures(pred_trainsmote))

    olr_corrected <- rbind(olr_corrected,
                           performanceMeasures(pred_train, 
                                               mean(datatrain$y)))
    
    bestth_rs <- threshold(pred_train$prediction, 
                           datatest$y)
    olr_bth_sp <- rbind(olr_bth_sp, 
                        performanceMeasures(pred_train, 
                                            bestth_rs))
    
    bestth_f1 <- threshold_f1(pred_train$prediction, 
                              datatest$y)
    olr_bth_f1 <- rbind(olr_bth_f1, 
                        performanceMeasures(pred_train, 
                                            bestth_f1))
    
    bestth_APR_curve <- threshold_APR_curve(pred_train$prediction, 
                                            datatest$y)
    olr_bth_APR_curve <- rbind(olr_bth_APR_curve, 
                               performanceMeasures(pred_train, bestth_APR_curve))
    
  }
  
  olr[is.na(olr)] <- 0
  out <- list(olr=olr, 
              ros=ros, 
              rus=rus, 
              smt=smt, 
              olr_corrected=olr_corrected, 
              olr_bth_sp=olr_bth_sp,
              olr_bth_f1=olr_bth_f1,
              olr_bth_APR_curve = olr_bth_APR_curve, 
              mean_y = meanY)
  return(out)
}

startSeed0 <- 111111

coef<-c(0.005,0.09,-0.6,0.002,0.75,0.69, .001, -.56, -.00087)
coef1 <- coef + c(1, 0, 0, 0,   0, 0, 0, 0, 0)
coef2 <- coef + c(0, 0, 0, 1.5, 0, 0, 0, 0, 0)
coef3 <- coef + c(0, 0, 0, 0,   0, 0, 0, 0, 2)
coef4 <- coef + c(1, 0, 0, 1.5, 0, 0, 0, 0, 0)
coef5 <- coef + c(1, 0, 0, 0,   0, 0, 0, 0, 2)
coef6 <- coef + c(0, 0, 0, 1.5, 0, 0, 0, 0, 2)
coef7 <- coef + c(1, 0, 0, 1.5, 0, 0, 0, 0, 2)

listcoef <- out <- list()
listcoef[[1]] <- coef
listcoef[[2]] <- coef1
listcoef[[3]] <- coef2
listcoef[[4]] <- coef3
listcoef[[5]] <- coef4
listcoef[[6]] <- coef5
listcoef[[7]] <- coef6
listcoef[[8]] <- coef7


if(use_cluster){
# Want to know the time
start <- proc.time()
# Make clusters to parallel
cl <- parallel::makeForkCluster(8)
doParallel::registerDoParallel(cl)

tmp <- list()

tmp <- foreach(prev= prevalence) %:% foreach(i=1:length(listcoef), .packages = c('ROCR', 'smotefamily', 'caret', 
                                                                                 'brglm2', 'RSBID', 'PRROC', 'pROC', 
                                                                                 'dplyr')) %dopar% {genSim(R, n, prev, listcoef[[i]], startSeed0)}

# Stop clusters.
parallel::stopCluster(cl)
print(proc.time() - start)
}

if(!use_cluster){
  start <- proc.time()
  
  tmp <- list()
  
  for(prev in prevalence){
    p <- paste('p-',prev, sep = '')
    print(p)
    tmp[[p]] <- list()
    
    for(i in 1:length(listcoef)){
      ii <- paste('S', i, sep = '')
      
      tmp[[p]][[ii]] <- genSim(R, n, prev, listcoef[[i]], startSeed0)
    }
  }
}

# May 2025; also save the Rdata.
save(tmp, file = paste('result-', Sys.Date(), '.Rdata', sep = ''))
