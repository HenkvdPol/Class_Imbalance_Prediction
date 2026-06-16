################################################################################
# Load the dataframe
################################################################################
Res <- paste0('./Setting/Result/', list.files('./Setting/Result/'))


sim_names <- lapply(Res, function(x){strsplit(x, '-', fixed = TRUE)[[1]][c(2,3,6)]})
sim_names <- matrix(unlist(sim_names), ncol = 3, byrow = TRUE)
colnames(sim_names) <- c('sample size', 'dim', 'sim number')


################################################################################
# Necessary Libraries
################################################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(patchwork)
library(ggrastr) 

################################################################################
# Colour coding
################################################################################
PRIME_ROSE_colors <- c(
  "#001F5F", # navy blue
  "#3A0D73", # dark purple
  '#702572', # PRIME-ROSE blue
  '#BD1A8D', # PRIME-ROSE purple
  "#D1554F", # warm red
  "#E17A39", # burnt orange
  "#F4991E", # orange
  '#F37021'  # PRIME-ROSE orange
)
PRIME_ROSE_gradient <- colorRampPalette(PRIME_ROSE_colors)

################################################################################
# Pre-processing
################################################################################

pre_processing <- function(df){
  # Some outcome values are NA, running 
  table(df$Measures[is.na(df$Value)])
  # shows this is only the case for the precision and recall. In these situations
  # we have a division by 0 error since no positive prediction is made. This 
  # means we are actually in a 0/0-situation and precision, F1-score is 0.
  df$Value[is.na(df$Value)] <- 0
  
  # Rename some model or measure to align the manuscript
  df$Model <- factor(df$Model, levels = c(levels(df$Model), 'Standard 0.5-threshold'))
  df$Model[which(df$Model == 'Original')] <- 'Standard 0.5-threshold'
  
  df$Measures <- factor(df$Measures, levels = c(levels(df$Measures), 'F1-score'))
  df$Measures[which(df$Measures == 'F1')] <- 'F1-score'
  
  # We factor the outcome of some variables to allow good ggploting. Also the order
  # of plotting
  df$Measures <- factor(df$Measures, levels = c('F1-score', 'Precision', 'Recall', 'Specificity',
                                                   'Accuracy', 'FNR', 'FPR', 'F1', 'AUC', 'APR'))
  df$Model <- factor(df$Model, levels = c('APR-threshold',
                                          'AUC-threshold',
                                          'F1-threshold',
                                          'Prevalence-threshold',
                                          'Standard 0.5-threshold',
                                          'RUS',
                                          'ROS',
                                          'SMOTE', 'Original'))
  df$Prevalence <- as.factor(df$Prevalence)
  df$Dimension <- as.factor(df$Dimension)
  df$delta <- as.factor(df$delta)
  df$n <- as.factor(df$n)
  df$Scenario <- as.factor(df$Scenario)
  
  return(df)
}

retrieve_dataframes <- function(measures_of_interest = c("Accuracy", "Precision", "Recall", "Specificity", 
                                                         "FNR", "FPR", "F1", "AUC", "APR", "F1-score"), 
                                models_of_interest = c("Original", "ROS", "RUS", "SMOTE", 
                                                       "Prevalence-threshold", "AUC-threshold", "F1-threshold", "APR-threshold",
                                                       "Standard 0.5-threshold", "F1"), 
                                dimension_of_interest = c(1,2,4,8,16), 
                                sample_size_of_interest = c(3000, 10000),
                                prevalence_of_interest = c(0.025, 0.05, 0.1, 0.2),
                                scenarios_of_interest = c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 2)){
  df <- NULL
  
  index_sample_size <- which(sim_names[, 'sample size'] %in% sample_size_of_interest)
  index_dimension <- which(sim_names[, 'dim'] %in% dimension_of_interest)
  
  index <- intersect(index_sample_size, index_dimension)
  
  for(i in index){
    load(Res[i])
    df <- rbind(df, Result)
  }
  
  df <- pre_processing(df)
  
  df_index <- which(df$Measures %in% measures_of_interest &
                      df$Model %in% models_of_interest &
                      df$Prevalence %in% prevalence_of_interest & 
                      df$Scenario %in% scenarios_of_interest)
  
  df <- df[df_index, ]
  
  return(df)
}

################################################################################
# Figure A
################################################################################

measures_of_interest <- c('F1-score', 'Precision', 'Recall', 'Specificity')
models_of_interest <- c('AUC-threshold', 'F1-threshold',
                       'Standard 0.5-threshold', 'Prevalence-threshold',
                       'RUS', 'ROS', 'SMOTE')

figure_A <- function(){
# In the first figure of the result, we limit ourselves to the 1-dimensional
# case. We now do n = 3000
measures_of_interest <- c('F1-score', 'Precision', 'Recall', 'Specificity')
models_of_interest <- c('AUC-threshold', 'F1-threshold',
                       'Standard 0.5-threshold', 'Prevalence-threshold',
                       'RUS', 'ROS', 'SMOTE')

print('Figure A')
for(d in c(1,2,4,8,16)){
  print(d)
  for(s in c(3000, 10000)){
    print(s)
    df_figure_a <- retrieve_dataframes(measures_of_interest, models_of_interest,
                                       dimension_of_interest = d,
                                       sample_size_of_interest = s)
    
    

    
    plot_figure_a <- ggplot(data = df_figure_a,
                            aes(x = Scenario, y = Value, group = Prevalence, color = Prevalence)) +
      stat_summary(data = df_figure_a, 
                   fun = 'mean', geom = 'line', linetype = 'dashed') +
      stat_summary(data = df_figure_a,
                   fun = 'mean', geom = 'point') +
      labs(x = expression("Value of "~ beta[1]),
           y = "Mean performance measure") + #, title = paste('Number of variables:', unique(df_new_figure_2$Dimension))) +
      scale_x_discrete(limits = rev,
                       labels = c('2', '1', '0'),
                       breaks = c('11', '6', '1')) +
      facet_grid(Measures~Model) +  
      scale_color_manual(values = c('0.025' = PRIME_ROSE_colors[1],
                                    '0.05' = PRIME_ROSE_colors[3],
                                    '0.1' = PRIME_ROSE_colors[6],
                                    '0.2' = PRIME_ROSE_colors[8]),
                         labels = c('0.025' = '2.5%',
                                    '0.05'= '5%',
                                    '0.1' = '10%',
                                    '0.2' = '20%')) + 
      theme(legend.position = 'top')
    
    ggsave(plot_figure_a, file = paste0('Result/FigureA-', s, '-', d, '.pdf'), width = 12.5, height = 14)
    rm(df_figure_a, plot_figure_a)
  }
}
}
################################################################################
# Figure B
################################################################################


figure_B <- function(N = c(3000, 10000)){
print('Figure B')
for(n in N){
  print(n)
  tmp <- retrieve_dataframes(measures_of_interest, models_of_interest,
                             sample_size_of_interest = n) 
  for(p in c(0.025, 0.05, 0.1, 0.2)){
    print(p)
    df_figure_b <- tmp[tmp$Prevalence == p, ]
    
    # Have the same order as in figure A
    df_figure_b$Measures <- factor(df_figure_b$Measures,
                                   levels = c('F1-score', 'Precision',
                                              'Recall', 'Specificity'))
    
    df_figure_b$`Covariate dimension` <- as.factor(df_figure_b$Dimension)
    
    plot_figure_b <- ggplot(data = df_figure_b,
                            aes(x = Scenario, y = Value, 
                                group = `Covariate dimension`, 
                                color = `Covariate dimension`)) +
      
      stat_summary(data = df_figure_b, 
                   fun = 'mean', geom = 'line', linetype = 'dashed') +
      stat_summary(data = df_figure_b,
                   fun = 'mean', geom = 'point') +
      
      labs(x = expression("Value of "~ beta[1]),
           y = "Mean performance measure") +
      scale_x_discrete(limits = rev,
                       labels = c('2', '1', '0'),
                       breaks = c('11', '6', '1')) +
      
      facet_grid(Measures~Model) +  
      scale_color_manual(values = c('1' = PRIME_ROSE_colors[1],
                                    '2' = PRIME_ROSE_colors[2],
                                    '4' = PRIME_ROSE_colors[4],
                                    '8' = PRIME_ROSE_colors[6],
                                    '16' = PRIME_ROSE_colors[8])) + 
      theme(legend.position = 'top')
    
    # Save the plot
    ggsave(plot_figure_b, file = paste0('Result/FigureB-',n,'-',p*100,'.pdf'), width = 12.5, height = 14)
    
    rm(df_figure_b, plot_figure_b)
  }
  rm(tmp)
}
}
################################################################################
# Figure C
################################################################################

figure_C <- function(){
print('Figure C')

# This figure is plotting the ROS, RUS and SMOTE calibration plots..
for(d in c(1,2,4,8,16)){
  print(d)
  for(n in c(3000, 10000)){
    print(n)
    tmp <- retrieve_dataframes(measures_of_interest = c("Accuracy", "Precision", "Recall", "Specificity", "FNR", "FPR", "F1", "AUC", "APR", "F1-score"),
                               models_of_interest = c("Original", "ROS", "RUS", "SMOTE", "Prevalence-threshold", 'Standard 0.5-threshold',
                                                        "AUC-threshold", "F1-threshold", "APR-threshold", 'F1'),
                               dimension_of_interest = d,
                               sample_size_of_interest = n)
    print(dim(tmp))
    for(m in c('ROS', 'RUS', "SMOTE", "Standard 0.5-threshold")){
      print(m)
      df_figure_C <- tmp[tmp$Model == m, ]
      
      df_figure_C$Model <- factor(df_figure_C$Model,
                                  levels = c(levels(df_figure_C$Model),
                                             c('x', 'y' , 'Calibration slope',
                                               'Calibration intercept', 'APR', 'AUC')))
      
      # We can do the following as the dataframe is created using expand.grid.
      df_figure_C$Model[which(df_figure_C$Measures == 'Precision')] <- 'x'
      df_figure_C$Model[which(df_figure_C$Measures == 'Specificity')] <- 'y'
      df_figure_C$Model[which(df_figure_C$Measures == 'FPR')] <- 'Calibration intercept'
      df_figure_C$Model[which(df_figure_C$Measures == 'APR')] <- 'APR'
      df_figure_C$Model[which(df_figure_C$Measures == 'AUC')] <- 'AUC'
      df_figure_C$Model[which(df_figure_C$Measures == 'F1-score')] <- 'Calibration slope'
      
      
      # Now change the outcome to the appropriate setting
      df_figure_C$Calibration[which(df_figure_C$Calibration == 'slope')] <- df_figure_C$Calibration[which(df_figure_C$Calibration == 'slope') + 1]
      df_figure_C$Calibration[which(df_figure_C$Measures == 'AUC')] <- df_figure_C$Value[which(df_figure_C$Measures == 'AUC')]
      df_figure_C$Calibration[which(df_figure_C$Calibration == 'check')] <- df_figure_C$Value[which(df_figure_C$Calibration == 'check')]
      
      
      # Clean up the dataframe.
      df_figure_C <- df_figure_C[-which(df_figure_C$Model %in% c(m, 'x', 'y')), ]
      
      df_figure_C <- df_figure_C[which(df_figure_C$Scenario %in% c('11', '6', '1')), ]
      df_figure_C$Calibration <- as.numeric(df_figure_C$Calibration)
      levels(df_figure_C$Prevalence) <- c('2.5%', '5%', '10%', '20%')
      
      # Remove extreme-values (larger than 10)
      index_to_remove <- which((df_figure_C$Calibration)^2 >= 10^2)
      if(length(index_to_remove) > 0){
        print(paste(length(index_to_remove), 'items removed.'))
        df_figure_C <- df_figure_C[-index_to_remove, ]
      }
      
      print(which((df_figure_C$Calibration[index])^2 >= 10^2))
      

      
      hlines <- data.frame(Model = c("Calibration slope", "Calibration intercept"), 
                           yintercept = c(1, 0))
      
      plot_figure_c <-  ggplot(data = df_figure_C,
                               aes(x = Scenario, y = Calibration, group = Scenario)) +
        geom_hline(data = hlines, aes(yintercept = yintercept), 
                   color = PRIME_ROSE_colors[8], linewidth = 1) +
        geom_boxplot() + 
        labs(x = expression("Value of "~ beta[1]),
             y = "Outcome") +
        facet_grid(Model~Prevalence, scales = 'free_y') + 
        scale_x_discrete(limits = rev,
                         labels = c('2', '1', '0'),
                         breaks = c('11', '6', '1'))
      ggsave(plot_figure_c, file = paste0('Result/FigureC-',m,'-',d, '-', n, '.pdf'), width = 12.5, height = 14)
      
      rm(df_figure_C, plot_figure_c)
    }}}
  
}
################################################################################
# Figure D
################################################################################
figure_D <- function(){
models_of_interest <- c('AUC-threshold', 'F1-threshold',
                       'Prevalence-threshold')
print('Figure D')
# This one is for showing the box-plots of the thresholds.

for(n in c(3000, 10000)){
  print(n)
  for(d in c(1,2,4,8,16)){
    print(d)
    tmp <- retrieve_dataframes(measures_of_interest = c("Precision", "Recall", "Specificity", "F1","F1-score"),
                               models_of_interest = c("Prevalence-threshold","AUC-threshold", "F1-threshold",'F1'),
                               dimension_of_interest = d,
                               sample_size_of_interest = n)
    
    for(p in c(0.025,0.05,0.1,0.2)){
      print(p)
      df_figure_d <- tmp[tmp$Prevalence == p, ]

plot_figure_d <- ggplot(data = df_figure_d,
                        aes(x = Scenario, y = Threshold)) +
  stat_summary(data = df_figure_d, fun = 'mean', geom = 'point') +
  stat_summary(
    data = df_figure_d,
    fun.data = function(x) data.frame(
      ymin   = quantile(x, 0.25),
      lower  = quantile(x, 0.25),
      middle = median(x),
      upper  = quantile(x, 0.75),
      ymax   = quantile(x, 0.75)
    ),
    geom = "boxplot",
    width = 0.5
  ) +
  labs(x = expression("Value of "~ beta[1]),
       y = "Threshold") +
  scale_x_discrete(limits = rev,
                   labels = c('2', '1', '0'),
                   breaks = c('11', '6', '1')) +
  facet_grid(Measures~Model) +  
  theme(legend.position = 'top')

ggsave(plot_figure_d, file = paste0('Result/FigureD-threshold-',n,'-', p*100, '-', d, '.pdf'), width = 12.5, height = 14)

rm(df_figure_d, plot_figure_d)
    }
  }
}
}
################################################################################
# Figure E
################################################################################
figure_E <- function(){
print('Figure E')
# This plots all the box-plots for all the measures
for(n in c(3000, 10000)){
  print(n)
  for(d in c(1,2,4,8,16)){
    print(d)
    tmp <- retrieve_dataframes(measures_of_interest = c("Precision", "Recall", "Specificity", "F1","F1-score"),
                               models_of_interest = c("ROS", "RUS", "SMOTE", "Prevalence-threshold", 
                                                      "AUC-threshold", "F1-threshold", 
                                                      'Standard 0.5-threshold'), 
                               dimension_of_interest = d,
                               sample_size_of_interest = n)
    
    for(p in c(0.025,0.05,0.1,0.2)){
      print(p)
      df_figure_e <- tmp[tmp$Prevalence == p, ]
      
      plot_figure_e <- ggplot(data = df_figure_e,
                              aes(x = Scenario, y = Value)) +
        stat_summary(data = df_figure_e, fun = 'mean', geom = 'point') +
        stat_summary(
          data = df_figure_e,
          fun.data = function(x) data.frame(
            ymin   = quantile(x, 0.25),
            lower  = quantile(x, 0.25),
            middle = median(x),
            upper  = quantile(x, 0.75),
            ymax   = quantile(x, 0.75)
          ),
          geom = "boxplot",
          width = 0.5
        ) +
        labs(x = expression("Value of "~ beta[1]),
             y = "Value of performance measure") +
        scale_x_discrete(limits = rev,
                         labels = c('2', '1', '0'),
                         breaks = c('11', '6', '1')) +
        facet_grid(Measures~Model) +  
        theme(legend.position = 'top')
      
      ggsave(plot_figure_e, file = paste0('Result/FigureE-',n,'-', p*100, '-', d, '.pdf'), width = 12.5, height = 14)
      
      rm(df_figure_e, plot_figure_e)
    }
  }
}
}

print('END')









