# What we are going to do here is, plot all figures and safe them as rplots. Then we are going
# to make one figure of multiple plots. This saves a lot of memory.

################################################################################
# Load the dataframes
################################################################################

Res <- paste0('./Setting/Result/', list.files('./Setting/Result/'))

################################################################################
# Necessary Libraries
################################################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(patchwork)
library(ggrastr) 
library(data.table)
library(patchwork)
library(grid)
library(gridExtra)

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

pre_processing <- function(df, model_to_use = c("Standard 0.5-threshold","ROS","RUS","SMOTE")){

# Rename some model or measure to align the manuscript
df$Model <- factor(df$Model, levels = c(levels(df$Model), 'Standard 0.5-threshold'))
df$Model[which(df$Model == 'Original')] <- 'Standard 0.5-threshold'

df$Measures <- factor(df$Measures, levels = c(levels(df$Measures), 'F1-score'))
df$Measures[which(df$Measures == 'F1')] <- 'F1-score'

# We now change the dataframe such that we only have the calibration
# plots including the different models, scenarios, prevalences, etc.

df$Measures <- factor(df$Measures,
                   levels = c(levels(df$Measures), 
                              c('x', 'y' ,
                                'Calibration slope','Calibration intercept')))

# We can do the following as the dataframe is created using expand.grid.
df$Measures[which(df$Measures == 'Precision')] <- 'x'
df$Measures[which(df$Measures == 'Specificity')] <- 'y'
df$Measures[which(df$Measures == 'FPR')] <- 'Calibration intercept'
df$Measures[which(df$Measures == 'AUC')] <- 'Calibration slope'

# Now change the outcome to the appropriate setting
df$Calibration[which(df$Calibration == 'slope')] <- df$Calibration[which(df$Calibration == 'slope') + 1]
df$Calibration[which(df$Measures == 'AUC')] <- df$Value[which(df$Measures == 'AUC')]
df$Calibration[which(df$Calibration == 'check')] <- df$Value[which(df$Calibration == 'check')]

# Clean up the dataframe.
df <- df[which(df$Measures %in% c('x', 'y')),]

# Only use the scenarios that is reasonable to plot.
# tadf <- df[which(df$Scenario %in% c('11', '6', '1')), ]
levels(df$Prevalence) <- c('2.5%', '5%', '10%', '20%')

# then we remove the calibration rows from the different thresholds.
df <- df[!is.na(df$Calibration),]

# And only keep the variables we need.
df <- df[, c('Measures', 'Model', 'Calibration', 'Prevalence', 
             'Dimension', 'n', 'Scenario', 'Simulation')]

# Maybe we want to use only one model.
df <- df[which(df$Model %in% model_to_use), ]

return(df)
}

################################################################################
# Expand dataframe.
################################################################################

make_plot <- function(long_dt, mean_dt, show_x = FALSE, show_y = FALSE) {
  p <- ggplot() +
    rasterise(geom_line(data = long_dt, aes(x = x, y = y, group = trace),
                        colour = "grey70", alpha = 0.3, linewidth = 0.2),
              dpi = 300) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),
                 colour = "red", linetype = "dashed", linewidth = 0.8) +
    geom_line(data = mean_dt, aes(x = x, y = y),
              colour = "black", linewidth = 1.2) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    scale_x_continuous(breaks = rep(0:10)/10, expand = expansion(mult = c(0.05, 0.05))) +
    scale_y_continuous(breaks = rep(0:10)/10, expand = expansion(mult = c(0.05, 0.05))) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          plot.margin = margin(10, 10, 10, 10))
  
  if(!show_x && !show_y){
    p <- p + theme(axis.text.x = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.ticks.x = element_blank())
  }
  
  if (show_x) {
    p <- p + theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.line.x = element_line(colour = "black", linewidth = 1))
  }
  if (show_y) {
    p <- p + theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.line.y = element_line(colour = "black", linewidth = 1))
  }
  
  if(show_x && show_y){
    p <- ggplot() +
      rasterise(geom_line(data = long_dt, aes(x = x, y = y, group = trace),
                          colour = "grey70", alpha = 0.3, linewidth = 0.2),
                dpi = 300) +
      geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),
                   colour = "red", linetype = "dashed", linewidth = 0.8) +
      geom_line(data = mean_dt, aes(x = x, y = y),
                colour = "black", linewidth = 1.2) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      scale_x_continuous(breaks = rep(0:10)/10, expand = expansion(mult = c(0.05, 0.05))) +
      scale_y_continuous(breaks = rep(0:10)/10, expand = expansion(mult = c(0.05, 0.05))) +
      theme_minimal() +
      theme(axis.title = element_blank(),
            plot.margin = margin(10, 10, 10, 10),
            axis.line.y = element_line(colour = "black", linewidth = 1),
            axis.line.x = element_line(colour = "black", linewidth = 1))
  }
  
  return(p)
}


manual_expand_to_plot <- function(df, n, show_x, show_y){
  
  n_traces <- length(df$Calibration) / 2
  
  # Split all strings in one shot
  parsed <- lapply(df$Calibration, function(s) as.numeric(strsplit(s, ",", fixed = TRUE)[[1]]))
  
  # Odd indices are x, even indices are y
  x_list <- parsed[seq(1, length(parsed), by = 2)]
  y_list <- parsed[seq(2, length(parsed), by = 2)]
  
  lengths_per_trace <- lengths(x_list)
  
  long_dt <- data.table(
    trace = rep(seq_len(n_traces), lengths_per_trace),
    x = unlist(x_list, use.names = FALSE),
    y = unlist(y_list, use.names = FALSE))
  
  common_x <- seq(min(long_dt$x), max(long_dt$x), length.out = ifelse(n = '3000', 
                                                                      1000,
                                                                      3333))
  
  interp_dt <- long_dt[, .(x = common_x,
                           y = approx(x, y, xout = common_x, rule = 2)$y), 
                       by = trace]
  
  mean_dt <- interp_dt[, .(y = mean(y)), by = x]
  
  p <- make_plot(long_dt, mean_dt, show_x, show_y)
  
  return(p)
}

# Using this we can now create all the total images.
strip_label <- function(text, pos){
  if(pos == 'horizontal'){
    p <- ggplot() +
      annotate("text", x = 1/2, y = 0.5, label = text,
               size = 4) +
      ylim(c(0,1)) + xlim(c(0,1)) +
      theme_void() +
      theme(plot.margin = margin(1/2, 1/2, 1/2, 1/2))
  }
  if(pos == 'vertical'){
    p <- ggplot() +
      annotate("text", x = 0.01, y = 0.5, label = text,
               size = 4) +
      ylim(c(0.45, 0.55)) + xlim(c(0, 0.15)) +
      theme_void() +
      theme(plot.margin = margin(1, 1, 1, 1))
  }
  
  return(p)
}

################################################################################
# Calibration plots
################################################################################

# For each dataframe, we retrieve the name, this indicates the sample size and
# dimension (in the name) and scenario and prevalence (in the number)
sim_numbers <- c(1,2,3,4,21,22,23,24,41,42,43,44)
sim_numbers <- paste0(sim_numbers, '.Rdata')

sim_names <- lapply(Res, function(x){strsplit(x, '-', fixed = TRUE)[[1]][c(2,3,6)]})
sim_names <- matrix(unlist(sim_names), ncol = 3, byrow = TRUE)
colnames(sim_names) <- c('sample size', 'dim', 'sim number')


calibration_data <- function(){
for(n in c("3000", "10000")){
#n = 10000
  print(n)
  for(d in c("1", "2", "4", "8", "16")){
  #d = 2
    print(d)
    index_sample_size <- which(sim_names[, 'sample size'] == n)
    index_dimension <- which(sim_names[, 'dim'] == d)
    
    index <- intersect(index_sample_size, index_dimension)
    
    for(s in sim_numbers){
      index_sim <- which(sim_names[index, 'sim number'] == s)
      
      # We now have the indexes of the result files.
      df <- NULL
      for(i in index_sim){
        load(Res[i])
        df <- rbind(df, Result)
      }
      
      for(m in c("Standard 0.5-threshold","ROS","RUS","SMOTE")){
      #m = 'SMOTE'
        df_plot <- pre_processing(df, model_to_use = m)
        
        print(m)
        
        show_y <- ifelse(s %in% paste0(c(1,2,3,4), '.Rdata'), TRUE, FALSE)
        show_x <- ifelse(s %in% paste0(c(4,24,44), '.Rdata'), TRUE, FALSE)
        
        p <- manual_expand_to_plot(df_plot, n, show_x, show_y)
        
        saveRDS(p, file = paste('Rplots/plot', n, d, gsub('.Rdata', '', s), m, '.rds', sep = '-'))
        
        rm(p, df_plot, show_x, show_y)
      }
      
      rm(df)
    }
  }
}
}

# Column strips (above the top row)
col1_strip <- strip_label(expression(beta[1] == 0), pos = 'horizontal')
col2_strip <- strip_label(expression(beta[1] == 1), pos = 'horizontal')
col3_strip <- strip_label(expression(beta[1] == 2), pos = 'horizontal')

# Row strips (to the right of each row), rotated 270°
row1_strip <- strip_label("2.5%", pos = 'vertical')
row2_strip <- strip_label("5%", pos = 'vertical')
row3_strip <- strip_label("10%", pos = 'vertical')
row4_strip <- strip_label("20%", pos = 'vertical')

# Empty spacer for corner cells of the strip header row
spacer <- plot_spacer()

# load all the individual plots.
calibration_plot <- function(){
for(n in c("3000", "10000")){
  print(n)
  for(d in c("1", "2", "4", "8", "16")){
    print(d)
    for(m in c("Standard 0.5-threshold","ROS","RUS","SMOTE")){
      print(m)
plot11 <- readRDS(file = paste('Rplots/plot', n, d, 1, m, '.rds', sep = '-'))
plot12 <- readRDS(file = paste('Rplots/plot', n, d, 21, m, '.rds', sep = '-'))
plot13 <- readRDS(file = paste('Rplots/plot', n, d, 41, m, '.rds', sep = '-'))
plot21 <- readRDS(file = paste('Rplots/plot', n, d, 2, m, '.rds', sep = '-'))
plot22 <- readRDS(file = paste('Rplots/plot', n, d, 22, m, '.rds', sep = '-'))
plot23 <- readRDS(file = paste('Rplots/plot', n, d, 42, m, '.rds', sep = '-'))
plot31 <- readRDS(file = paste('Rplots/plot', n, d, 3, m, '.rds', sep = '-'))
plot32 <- readRDS(file = paste('Rplots/plot', n, d, 23, m, '.rds', sep = '-'))
plot33 <- readRDS(file = paste('Rplots/plot', n, d, 43, m, '.rds', sep = '-'))
plot41 <- readRDS(file = paste('Rplots/plot', n, d, 4, m, '.rds', sep = '-'))
plot42 <- readRDS(file = paste('Rplots/plot', n, d, 24, m, '.rds', sep = '-'))
plot43 <- readRDS(file = paste('Rplots/plot', n, d, 44, m, '.rds', sep = '-'))


combined <- (col1_strip | col2_strip | col3_strip | spacer) /
  (plot11    | plot12        | plot13        | row1_strip) /
  (plot21    | plot22        | plot23        | row2_strip) /
  (plot31    | plot32        | plot33        | row3_strip) /
  (plot41    | plot42        | plot43        | row4_strip) +
  plot_layout(heights = c(0.1, 1, 1, 1, 1),
              widths  = c(1, 1, 1, 0.05))

final_wrapped <- grid.arrange(
  patchworkGrob(combined),
  left = textGrob("Observed proportion", rot = 90,
                  gp = gpar(fontsize = 13)),
  bottom = textGrob("Estimated probabilities", hjust = 1,
                    gp = gpar(fontsize = 13))
)

final_wrapped

filename <- paste('Result/Calibration-plot', n, d, m, '.pdf', sep = '-')

ggsave(filename, final_wrapped, width = 11, height = 9)
system(paste0("pdfcrop ", shQuote(filename), " ", shQuote(filename)))

rm(final_wrapped, combined, plot11, plot12, plot13,
   plot21, plot22, plot23, plot31, plot32, plot33)
    }}}
}


