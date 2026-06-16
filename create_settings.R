# Here we create all the settings.

source('Setting/setting.R')

seeds <- c(0, 1234, 12345, 123456)

parameters <- expand.grid('n' = n,
                          'dimension' = dimension,
                          'delta' = delta,
                          'seed' = seeds)

dir.create('Setting', showWarnings = FALSE, recursive = TRUE)
dir.create('Setting/Result', showWarnings = FALSE, recursive = TRUE)

total_text <- NULL
for(i in nrow(parameters):1){
  n <- as.integer(parameters$n[i])
  p <- as.integer(parameters$dimension[i])
  d <- as.double(parameters$delta[i])
  s <- as.integer(parameters$seed[i])
  
  print(total_text)
  sbatch_text <- paste0( n, ' ', p, ' ', d, ' ', s, ' ',use_cluster, ' ', R)
  total_text <- c(total_text, sbatch_text)
  
}
writeLines(total_text, file.path("Setting", paste0("parameters.txt")), sep = "\n", useBytes = TRUE)
