# Description

# Author: Timothy Keyes
# Version: 2020-02-05

# Libraries
library(flowCore)
library(tidyverse)

# source aml utils
source(here::here("scripts", "setup", "aml_utils.R"))

# source tidyTOF
source_tidyTOF(file.path("~", "GitHub", "tidyTOF"))

# Parameters

input_path <-
  file.path(
    "~", 
    "Box", 
    "Tim", 
    "Lab", 
    "Data", 
    "DDPR_Data" 
  )
out_data <- here::here("data")
marker_path <- here::here("docs", "ALL_panel.csv")

#===============================================================================

#read in names of markers in the dataset 
marker_names <-
  marker_path %>%
  read_csv() %>% 
  mutate(Metal = str_replace_all(Metal, "[()]", ""))
  

my_data <- 
  input_path %>%
  list.files(path = ., full.names = TRUE) %>% 
  str_split(pattern = "_", simplify = TRUE) %>% 
  as_tibble() %>% 
  transmute(
    file_name = list.files(path = input_path, full.names = TRUE),
    patient = str_split(V2, pattern = "/") %>% 
      map_chr(.f = last),
    stimulation = str_replace(V3, ".fcs", ""), 
    data = 
      map(
        file_name, 
        ~ 
          read.FCS(
            filename = ., 
            transformation = FALSE, 
            truncate_max_range = FALSE
          ) %>% 
          flowCore::exprs() %>% 
          as_tibble()
      )
  )

col_names <- 
  map(my_data$data, colnames) %>% 
  unlist() %>% 
  str_replace_all(pattern = "[()]", replacement = "") %>% 
  table() %>% 
  enframe()


#may have to change this data structure to use regular expressions
lookup_table <-
  setNames(object = marker_names$Metal, nm = marker_names$protein)

tof_rename <- function(data, lookup_table) { 
  colnames(data) <- str_replace_all(colnames(data), "[()]", "")
  my_lookups <- (lookup_table[which(lookup_table %in% colnames(data))])
  

  data %>% 
    select(which(colnames(.) %in% lookup_table)) %>% 
    rename(!!! my_lookups)
}

my_data <- 
  my_data %>% 
  mutate(num_cols = map_dbl(data, ncol)) %>% 
  dplyr::filter(num_cols > 1) %>% 
  mutate(
    data = map(data, tof_rename, lookup_table)
  ) %>% 
  select(-num_cols)

my_data <- 
  my_data %>% 
  group_by(file_name, patient, stimulation) %>% 
  unnest(cols = data) %>% 
  ungroup() %>% 
  tof_preprocess()

write_rds(
  x = my_data, path = file.path(out_data, "all_data_p.rds"), compress =  "gz"
)

sampled_data <- 
  my_data %>% 
  group_by(patient, stimulation) %>% 
  sample_n(size = 2000, replace = TRUE) %>% 
  ungroup()

write_rds(
  x = sampled_data, path = file.path(out_data, "sampled_all_data_p.rds"), compress =  "gz"
)
