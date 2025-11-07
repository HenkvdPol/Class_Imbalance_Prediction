# This file will create the plots for the report. For this we need the
# results for N=1000, and N=10000.
load('result-2025-09-19.Rdata')
result103 <- tmp
load('result-2025-09-24.Rdata')
result104 <- tmp
rm(tmp)

library(ggplot2)
library(dplyr)

################################################################################
#                    Preparing the dataframe
################################################################################

# First I need to merge everything in one list. Moreover, we want to have 
# appropriate names for the plots.
Models <- c('Standard 0.5-threshold', 'ROS', 'RUS', 'SMOTE', 'Prevalence-threshold',
            'AUC-threshold', 'F1-threshold', 'APR-threshold', 'PY')
Prevalence <- c(0.025, 0.05, 0.10, 0.20)
Scenarios <- (1:8)
Measures <- c('Accuracy', 'Precision', 'Recall', 'Specificity', 'FNR',
              'FPR', 'F1', 'AUC', 'APR')
Resamples <- 1:1000
sample_size <- c(1000, 10000)

result <- expand.grid(Resample = Resamples,
                      Measures = Measures,
                      Models = Models,
                      Scenarios = Scenarios,
                      Prevalence = Prevalence,
                      Sample_size = sample_size)

result$value <- NA

# the result(sample_size) dataframes are list of lists. In the for-loop
# they were not given appropriate names.
names(result104) <- Prevalence
names(result103) <- Prevalence

# Now start the for-loop to insert the results in one dataframe. We also
# name here the nested lists and dataframes. 
for(prev in names(result104)){
  names(result104[[prev]]) <- Scenarios
  names(result103[[prev]]) <- Scenarios
  
  for(scenario in names(result104[[prev]])){
    names(result103[[prev]][[scenario]]) <- Models
    names(result104[[prev]][[scenario]]) <- Models
    
    for(model in names(result104[[prev]][[scenario]][-9])){
      # column 9 shows the actual prevalence in the simulation. We omit it here
      colnames(result104[[prev]][[scenario]][[model]]) <- Measures
      colnames(result103[[prev]][[scenario]][[model]]) <- Measures
      
      
      for(measure in colnames(result104[[prev]][[scenario]][[model]])){
        print(paste(prev, scenario, model, measure))
        
        # Retrieve the values
        value103 <- result103[[prev]][[scenario]][[model]][, measure]
        value104 <- result104[[prev]][[scenario]][[model]][, measure]
        
        # Retrieve the indexes using which(). 
        # This is extremely slow as we have a very large list.
        index103 <- which(result$Prevalence == prev & result$Measures == measure & result$Models == model & result$Scenarios == scenario &  result$Sample_size == 1000)
        index104 <- which(result$Prevalence == prev & result$Measures == measure & result$Models == model & result$Scenarios == scenario & result$Sample_size == 10000)
        
        
        # Update the result values.
        result$value[index103] <- value103
        result$value[index104] <- value104
      }
    }
  }
}

# Some measures may be NA. In that case, it should be 0.
result$value[is.na(result$value)] <- 0

# The scenarios are now integer values. We need to change it to appropriate Si's.
result$Scenarios <- paste('S', result$Scenarios, sep = '')

# Finally, we won't show all the measures and models.
models_of_interest <- c('Standard 0.5-threshold', 'ROS', 'RUS', 'SMOTE', 'AUC-threshold',
                       'Prevalence-threshold', 'APR-threshold',
                       'F1-threshold')
measures_of_interest <- c('Precision', 'Recall', 'F1-score', 'AUC', 'APR')

# 17-10-2025; switch S4 and S5
S4 <- which(result$Scenarios == 'S4')
S5 <- which(result$Scenarios == 'S5')

result$Scenarios[S4] <- 'S5'
result$Scenarios[S5] <- 'S4'

# 29-10-2025; Last improvements to have consistent figures in manuscript.
result$Models <- factor(result$Models, levels = c('Standard 0.5-threshold', levels(result$Models)))
result$Measures <- factor(result$Measures, levels = c('F1-score', levels(result$Measures)))

result$Models[which(result$Models == 'Standard')] <- 'Standard 0.5-threshold'
result$Measures[which(result$Measures == 'F1')] <- 'F1-score'

result <- result[-which(result$Models == 'PY'),]

rm(S4, S5)
################################################################################
#                    Create Plots
################################################################################


# Colour gradients we will use to create the graphs.
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


Figure3 <- ggplot(data = result[which(result$Measures %in% measures_of_interest &
                                        result$Models %in% models_of_interest &
                                        result$Sample_size == 1000 &
                                        result$Scenarios == 'S1'), ], 
                  aes(x = Prevalence, y = value, group = Prevalence)) +
  geom_boxplot(color = PRIME_ROSE_colors[1]) + 
  ylab('Performance measure') +
  xlab('Prevalence') +
  facet_grid(Measures~Models, scales = 'free_y') +
  scale_x_continuous(breaks = c(0.025, 0.05, 0.1, 0.2), limits = c(0, 0.225),
                     labels = scales::percent_format(accuracy = 0.1) ) + 
  theme(legend.position = 'right', 
        axis.text.x = element_text(angle = -90, vjust = 0, hjust = 0, size = 8))
Figure3
ggsave(Figure3, file = 'Figure3.pdf', width = 12.5, height = 14)



Figure2 <- ggplot(data = result[which(result$Measures %in% measures_of_interest &
                                        result$Models %in% models_of_interest &
                                        result$Scenarios == 'S1'), ], 
                  aes(x = Prevalence, y = value, group = Sample_size)) +
  stat_summary(data = result[which(result$Measures %in% measures_of_interest & 
                                     result$Models %in% models_of_interest &
                                     result$Sample_size == 1000 &
                                     result$Scenarios == 'S1'),],
               aes(color = 'N = 1000'),
               fun = 'mean', geom = 'line', linetype = 'dashed') + #, color = PRIME_ROSE_colors[1]) +
  stat_summary(data = result[which(result$Measures %in% measures_of_interest & 
                                     result$Models %in% models_of_interest &
                                     result$Sample_size == 1000 &
                                     result$Scenarios == 'S1'),],
               aes(color = 'N = 1000'),
               fun = 'mean', geom = 'point') + #, color = PRIME_ROSE_colors[1]) +
  
  stat_summary(data = result[which(result$Measures %in% measures_of_interest & 
                                     result$Models %in% models_of_interest &
                                     result$Sample_size == 10000 &
                                     result$Scenarios == 'S1'),],
               aes(color = 'N = 10000'),
               fun = 'mean', geom = 'line', linetype = 'dashed',) + # color = PRIME_ROSE_colors[8]) +
  stat_summary(data = result[which(result$Measures %in% measures_of_interest & 
                                     result$Models %in% models_of_interest &
                                     result$Sample_size == 10000 &
                                     result$Scenarios == 'S1'),],
               aes(color = 'N = 10000'),
               fun = 'mean', geom = 'point') + #, color = PRIME_ROSE_colors[8]) +
  ylab('Mean performance measure') +
  xlab('Prevalence') +
  facet_grid(Measures~Models, scales = 'free_y') +
  
  scale_color_manual(name = "Sample size",
                     values = c("N = 1000" = "#001F5F", "N = 10000" = "#F37021")) +
  
  scale_x_continuous(breaks = c(0.025, 0.05, 0.1, 0.2), limits = c(0, 0.225),
                     labels = scales::percent_format(accuracy = 0.1) ) + 
  theme(legend.position = 'top', 
        axis.text.x = element_text(angle = -90, vjust = 0, hjust = 0, size = 8))
Figure2
ggsave(Figure2, file = 'Figure2.pdf', width = 12.5, height = 14)


FigureSM5 <- ggplot(data = result[which(result$Measures %in% measures_of_interest &
                                        result$Models %in% models_of_interest &
                                        result$Sample_size == 10000 &
                                        result$Scenarios == 'S1'), ], 
                  aes(x = Prevalence, y = value, group = Prevalence)) +
  geom_boxplot(color = PRIME_ROSE_colors[8]) + 
  ylab('Performance measure') +
  xlab('Prevalence') +
  facet_grid(Measures~Models, scales = 'free_y') +
  scale_x_continuous(breaks = c(0.025, 0.05, 0.1, 0.2), limits = c(0, 0.225),
                     labels = scales::percent_format(accuracy = 0.1) ) + 
  theme(legend.position = 'right', 
        axis.text.x = element_text(angle = -90, vjust = 0, hjust = 0, size = 8))
FigureSM5
ggsave(FigureSM5, file = 'FigureSM5.pdf', width = 12.5, height = 14)



FigureSM3 <- ggplot(data = result[which(result$Measures %in% measures_of_interest &
                                        result$Prevalence == 0.025), ], 
                  aes(x = Scenarios, y = value, group = Sample_size)) +
  stat_summary(data = result[which(result$Measures %in% measures_of_interest & 
                                     result$Sample_size == 1000 &
                                     result$Prevalence == 0.025),],
               aes(color = 'N = 1000'),
               fun = 'mean', geom = 'line', linetype = 'dashed') +
  stat_summary(data = result[which(result$Measures %in% measures_of_interest & 
                                     result$Sample_size == 1000 &
                                     result$Prevalence == 0.025),], 
               aes(color = 'N = 1000'),
               fun = 'mean', geom = 'point') +
  
  stat_summary(data = result[which(result$Measures %in% measures_of_interest & 
                                     result$Sample_size == 10000 &
                                     result$Prevalence == 0.025),], 
               aes(color = 'N = 10000'),
               fun = 'mean', geom = 'line', linetype = 'dashed') +
  stat_summary(data = result[which(result$Measures %in% measures_of_interest & 
                                     result$Sample_size == 10000 &
                                     result$Prevalence == 0.025),], 
               aes(color = 'N = 10000'),
               fun = 'mean', geom = 'point') +
  ylab('Mean Performance Measure') +
  xlab('Scenarios') +
  
  scale_color_manual(name = "Sample size",
                     values = c("N = 1000" = "#001F5F", "N = 10000" = "#F37021")) +
  
  facet_grid(Measures~Models, scales = 'free_y') +
  
  theme(legend.position = 'top', 
        axis.text.x = element_text(angle = -90, vjust = 0, hjust = 0, size = 8))
FigureSM3
ggsave(FigureSM3, file = 'FigureSM3.pdf', width = 12.5, height = 14)


FigureSM4 <- ggplot(data = result[which(result$Measures %in% measures_of_interest &
                                        result$Prevalence == 0.025 &
                                        result$Sample_size == 10000), ], 
                  aes(x = Scenarios, y = value, group = Scenarios)) +
  geom_boxplot(color = PRIME_ROSE_colors[8]) +
  ylab('Performance measure') +
  xlab('Scenarios') +
  facet_grid(Measures~Models, scales = 'free_y')
FigureSM4
ggsave(FigureSM4, file = 'FigureSM4.pdf', width = 12.5, height = 14)



FigureSM1 <- ggplot(data = result[which(result$Measures %in% measures_of_interest &
                                             result$Prevalence == 0.025 &
                                             result$Sample_size == 1000), ], 
                       aes(x = Scenarios, y = value, group = Models, color = Models, by = Models, fill = Models)) +
  stat_summary(fun = 'mean', geom = 'line', linetype = 'dashed', linewidth = .6, show.legend = FALSE, colour = PRIME_ROSE_colors[1]) +
  stat_summary(fun = 'mean', geom = 'point', show.legend = FALSE, size = 1.5, colour = PRIME_ROSE_colors[1]) +
  
  ylab('Mean performance measure') +
  xlab('Scenarios') + theme(legend.position = 'top', legend.text = element_text(size = 10), legend.title = element_text(size = 10), legend.key.size = unit(1.3, 'cm')) +
  facet_grid(Measures~Models, scales = 'free_y')
FigureSM1
ggsave(FigureSM1, file = 'FigureSM1.pdf', width = 18, height = 22)


FigureSM2 <- ggplot(data = result[which(result$Measures %in% measures_of_interest &
                                                  result$Models %in% models_of_interest &
                                                  result$Prevalence == 0.025 &
                                                  result$Sample_size == 1000), ], 
                            aes(x = Scenarios, y = value)) +
  geom_boxplot(color = PRIME_ROSE_colors[1]) + 
  ylab('Mean performance measure') +
  xlab('Scenarios') + theme(legend.position = 'top', legend.text = element_text(size = 10), legend.title = element_text(size = 10), legend.key.size = unit(1.3, 'cm')) +
  facet_grid(Measures~Models, scales = 'free_y')
FigureSM2
ggsave(FigureSM2, file = 'FigureSM2.pdf', width = 18, height = 22)


measures_of_interest <- c('Precision', 'Recall', 'F1-score')
scenario_models_of_interest <- c('Prevalence-threshold', 'Standard 0.5-threshold', "AUC-threshold" , "F1-threshold",  "APR-threshold")

Figure4 <- ggplot(data = result[which(result$Measures %in% measures_of_interest &
                                                  result$Models %in% scenario_models_of_interest &
                                                  result$Prevalence == 0.025 &
                                                  result$Sample_size == 1000), ], 
                            aes(x = Scenarios, y = value, group = Models, color = Models, by = Models, fill = Models)) +
  stat_summary(fun = 'mean', geom = 'line', linetype = 'dashed', linewidth = .6, show.legend = FALSE, colour = PRIME_ROSE_colors[1]) +
  stat_summary(fun = 'mean', geom = 'point', show.legend = FALSE, size = 1.5, colour = PRIME_ROSE_colors[1]) +
  ylab('Mean Performance Measure') +
  xlab('Scenarios') + theme(legend.position = 'top', legend.text = element_text(size = 10), legend.title = element_text(size = 10), legend.key.size = unit(1.3, 'cm')) +
  facet_grid(Measures~Models, scales = 'free_y')
Figure4
ggsave(Figure4, file = 'Figure4.pdf', width = 18, height = 22)

Figure5 <- ggplot(data = result[which(result$Measures %in% measures_of_interest &
                                                  result$Models %in% scenario_models_of_interest &
                                                  result$Prevalence == 0.025 &
                                                  result$Sample_size == 1000), ], 
                            aes(x = Scenarios, y = value)) +
  geom_boxplot(color = PRIME_ROSE_colors[1]) + 
  ylab('Mean performance measure') +
  xlab('Scenarios') + theme(legend.position = 'top', legend.text = element_text(size = 10), legend.title = element_text(size = 10), legend.key.size = unit(1.3, 'cm')) +
  facet_grid(Measures~Models, scales = 'free_y')
Figure5
ggsave(Figure5, file = 'Figure5.pdf', width = 18, height = 22)



measures_of_interest <- c('Precision', 'Recall', 'F1-score', 'AUC', 'APR')

FigureSM4 <- ggplot(data = result[which(result$Measures %in% measures_of_interest &
                                                  result$Models %in% models_of_interest &
                                                  result$Prevalence == 0.025 &
                                                  result$Sample_size == 10000), ], 
                            aes(x = Scenarios, y = value)) +
  geom_boxplot(color = PRIME_ROSE_colors[8]) + 
  ylab('Mean performance measure') +
  xlab('Scenarios') + theme(legend.position = 'top', legend.text = element_text(size = 10), legend.title = element_text(size = 10), legend.key.size = unit(1.3, 'cm')) +
  facet_grid(Measures~Models, scales = 'free_y')
FigureSM4
ggsave(FigureSM4, file = 'FigureSM4.pdf', width = 18, height = 22)



















