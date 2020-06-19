# aml_read_data.R 

### Description:
### A script that reads in the data for the aml project (from its raw form). 
### Input: Directory for the aml dataset (a folder full of .fcs files)
### Output: A tibble containing every cell from the aml project.

########################

# Libraries
library(flowCore)
library(tidyverse)

# source files
source(here::here("scripts", "setup", "aml_utils.R"))
source(file.path("~", "GitHub", "tidyTOF", "tof_read_fcs.R"))

# Parameters
set_global_variables(locale = "galaxia")

#set up directories
message("Setting Directory Information")
reader_plot_output <- here::here("figures", "processing")
if (!dir.exists(reader_plot_output)) dir.create(reader_plot_output, showWarnings = FALSE)

#read files
message("Reading in data...")

my_data <- 
  RAW_DATA_DIRECTORY %>% 
  tof_read_fcs(folder_path = .)

##########clean the data that was read in 
message("Cleaning data...")
my_data <- 
  my_data %>% 
  mutate(
    plate = 
      str_extract(file_names, pattern = "Plate[:digit:]+") %>% 
      str_to_lower(), 
    patient = 
      str_extract(file_names, pattern = "_([:alnum:]+)") %>% 
      str_sub(start = 2L), 
    stimulation = 
      str_extract(file_names, pattern = "[:space:][[:alpha:][:digit:]]+_") %>% 
      str_sub(start = 2L, end = -2L),
    condition = 
      str_extract(file_names, pattern = "(Dx)|(Rx)") %>%
      replace_na(replace = "healthy") %>% 
      str_to_lower()
  ) %>%
  select(
    -Time, 
    -Event_length, 
    -contains("Pd"), 
    -Center, 
    -Offset, 
    -Width, 
    -Residual, 
    -beadDist
  ) %>% 
  rename_with(
    .fn = 
      function(x) str_extract(x, pattern = "_.+") %>% 
      str_sub(start = 2L), 
    .cols = contains("_")
  ) %>% 
  rename(file_names = names)

message("Saving Files")

#save a sampled version of the giant dataframe for testing locally 
my_data %>% 
  write_rds(path = file.path(DATA_DIRECTORY, "raw_aml_tibble.rds"))

my_data %>% 
  group_by(file_names) %>% 
  slice_sample(prop = 0.1) %>% 
  write_rds(path = file.path(DATA_DIRECTORY, "raw_aml_tibble_sampled.rds"))


########provide some simple figures regarding the quality of the data

cell_counts <- 
  my_data %>% 
  count(patient, stimulation, condition)

message("Plotting quality report...")

cell_counts_plot <- 
  cell_counts %>% 
  mutate(
    patient = fct_reorder(patient, n), 
    condition = factor(condition, levels = c("healthy", "dx", "rx"))
  ) %>% 
  ggplot(aes(x = patient, y = n, color = stimulation, shape = condition)) + 
  geom_point(size = 3) + 
  scale_y_continuous(
    minor_breaks = NULL, 
    labels = scales::label_number(accuracy = 1, scale = 1/1e3, suffix = "K")
  ) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9)) + 
  labs(
    x = "Patient ID", 
    y = "Number of cells (thousands)", 
    color = "Stimulation ", 
    shape = "Condition"
  )

ggsave(
  filename = "cell_counts.pdf", 
  plot = cell_counts_plot, 
  device = "pdf", 
  path = READER_OUTPUT, 
  width = 7, 
  height = 6, 
  units = "in"
)

cell_counts_by_plate <- 
  my_data %>% 
  count(plate, patient, condition) %>% 
  mutate(plate = factor(plate, levels = str_c("plate", 1:16))) %>% 
  ggplot(aes(x = plate, y = n)) + 
  geom_point(size = 3, alpha = 0.6) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9)) + 
  labs(
    subtitle = "Each dot represents a unique sample tube",
    x = "Plate", 
    y = "Number of cells (thousands)"
  )

ggsave(
  filename = "cell_counts_by_plate.pdf", 
  plot = cell_counts_by_plate, 
  device = "pdf", 
  path = READER_OUTPUT, 
  width = 7, 
  height = 5, 
  units = "in"
)


