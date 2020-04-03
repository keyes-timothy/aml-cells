# Description: 
### This script will read in files from an aml data matrix and return a smaller 
### matrix that preserves similar features, but is much smaller in size. 

# Author: Timothy Keyes
# Version: 2020-03-31

# Libraries
library(tidyverse)

# Parameters
input_path <- here::here("data", "AML_matrix_clustered.rds")

output_path <- here::here("data", "AML_matrix_clustered_sampled.rds")

#===============================================================================

aml_data <- 
  input_path %>% 
  read_rds()

aml_sampled <- 
  aml_data %>% 
  group_by(patient, stimulation, condition, Mah.cluster) %>% 
  sample_n(size = 1000, replace = TRUE) %>% 
  distinct() %>% 
  filter(n() > 100) %>% 
  ungroup()

aml_sampled %>% 
  write_rds(path = output_path)
  
