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
  df$Measures <- as.factor(df$Measures)
  df$Model <- factor(df$Model, levels = c('Prevalence-threshold',
                                          'Standard 0.5-threshold', 
                                          'F1-threshold',
                                          'AUC-threshold',
                                          'APR-threshold',
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

retrieve_dataframes <- function(measures_of_interest = c("Accuracy", "Precision", "Recall", "Specificity", "FNR", "FPR", "F1", "AUC", "APR"), 
                                models_of_interest = c("Original", "ROS", "RUS", "SMOTE", "Prevalence-threshold", "AUC-threshold", "F1-threshold", "APR-threshold"), 
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
# Figure 2
################################################################################

figure_2 <- function(){
# In the first figure of the result, we limit ourselves to the 1-dimensional
# case.
measures_of_interest <- c('F1-score', 'Precision', 'Recall', 'Specificity')
#measures_of_interest <- c('Precision')
models_of_interest <- c('AUC-threshold', 'F1-threshold',
                       'Standard 0.5-threshold', 'Prevalence-threshold')


df_figure_2 <- retrieve_dataframes(measures_of_interest, models_of_interest,
                                       1, 10000)

# Have the same order as figure 1
df_figure_2$Measures <- factor(df_figure_2$Measures,
                               levels = c('F1-score', 'Precision',
                                          'Recall', 'Specificity'))

# Create a highlight for the figure 1 models
highlight <- data.frame(Model = c('Prevalence-threshold', 'Standard 0.5-threshold'))

plot_figure_2 <- ggplot(data = df_figure_2,
                        aes(x = Scenario, y = Value, group = Prevalence, color = Prevalence)) +
  geom_rect(data = highlight,
            aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "#DDEEFF", alpha = 0.5,
            inherit.aes = FALSE) +
  stat_summary(data = df_figure_2, 
               fun = 'mean', geom = 'line', linetype = 'dashed') +
  stat_summary(data = df_figure_2,
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
plot_figure_2

# Save the plot
ggsave(plot_figure_2, file = 'Result/Figure2.pdf', width = 12.5, height = 14)



# Delete the data not needed anymore.
rm(df_figure_2, plot_figure_2, index)
}
################################################################################
# Figure 3
################################################################################

figure_3 <- function(){
# In figure 3, we focus on the results when the dimension increases. We keep
# the prevalence at 5%.
  measures_of_interest <- c('F1-score', 'Precision', 'Recall', 'Specificity')
  #measures_of_interest <- c('Precision')
  models_of_interest <- c('AUC-threshold', 'F1-threshold',
                         'Standard 0.5-threshold', 'Prevalence-threshold')
  
df_figure_3 <- retrieve_dataframes(measures_of_interest, 
                                   models_of_interest,
                                   sample_size_of_interest =  10000,
                                   prevalence_of_interest = 0.05)


# Have the same order as figure 1
df_figure_3$Measures <- factor(df_figure_3$Measures,
                               levels = c('F1-score', 'Precision',
                                          'Recall', 'Specificity'))

# Create a highlight for the figure 1 models
highlight <- data.frame(Model = c('Prevalence-threshold', 'Standard 0.5-threshold'))

df_figure_3$`Covariate dimension` <- as.factor(df_figure_3$Dimension)

plot_figure_3 <- ggplot(data = df_figure_3,
                        aes(x = Scenario, y = Value, 
                            group = `Covariate dimension`, 
                            color = `Covariate dimension`)) +
  geom_rect(data = highlight,
            aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "#DDEEFF", alpha = 0.5,
            inherit.aes = FALSE) +
  
  
  stat_summary(data = df_figure_3, 
               fun = 'mean', geom = 'line', linetype = 'dashed') +
  stat_summary(data = df_figure_3,
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
plot_figure_3

# Save the plot
ggsave(plot_figure_3, file = 'Result/Figure3.pdf', width = 12.5, height = 14)

rm(plot_figure_3, df_figure_3, highlight)
}
################################################################################
# Figure 5
################################################################################

figure_5 <- function(){
df_figure_5 <- retrieve_dataframes(models_of_interest = 'Standard 0.5-threshold',
                                   sample_size_of_interest = 10000)

df_figure_5$Model <- factor(df_figure_5$Model,
                            levels = c(levels(df_figure_5$Model),
                                       c('x', 'y' , 'Calibration slope',
                                         'Calibration intercept', 'APR', 'AUC')))

# We can do the following as the dataframe is created using expand.grid.
df_figure_5$Model[which(df_figure_5$Measures == 'APR')] <- 'APR'
df_figure_5$Model[which(df_figure_5$Measures == 'AUC')] <- 'AUC'

df_figure_5 <- df_figure_5[which(df_figure_5$Measures %in% c('APR', 'AUC')),]

levels(df_figure_5$Prevalence) <- c('2.5%', '5%', '10%', '20%')

plot_figure_5 <- ggplot(data = df_figure_5,
                        aes(x = Scenario, y = Value, group = Dimension, color = Dimension)) +
  
  stat_summary(data = df_figure_5, 
               fun = 'mean', geom = 'line', linetype = 'dashed') +
  stat_summary(data = df_figure_5,
               fun = 'mean', geom = 'point') +
  labs(x = expression("Value of "~ beta[1]),
       y = "Mean performance measure") +
  ylim(c(0,1)) +
  scale_x_discrete(limits = rev,
                   labels = c('2', '1', '0'),
                   breaks = c('11', '6', '1')) +
  facet_grid(Model~Prevalence) +  
  scale_color_manual(values = c('1' = PRIME_ROSE_colors[1],
                                '2' = PRIME_ROSE_colors[2],
                                '4' = PRIME_ROSE_colors[4],
                                '8' = PRIME_ROSE_colors[6],
                                '16' = PRIME_ROSE_colors[8])) +  
  theme(legend.position = 'top')
plot_figure_5

ggsave(plot_figure_5, file = 'Result/Figure5.pdf', width = 12.5, height = 14 / 2)
}