###     aml_utils.R
# Description: 
### Sources script files in the aml-cells project repository for a variety of 
### helper functions. 

# Author: Timothy Keyes
# Version: 2020-06-11

#### call_libraries ####
# Description: A function to load and attach all libraries for an analysis. 
#
# Inputs: 
#     -libraries = a character vector of libraries that need to be loaded
#
# Outputs: None
#
# Side-effects: All libraries in `libraries` are loaded and attached to the session. 
call_libraries <- function(libraries) { 
  for (library in libraries) {
    # install from CRAN if not already present
    if (!require(library, character.only = TRUE)) {
      install.packages(library, dep = TRUE, quiet = TRUE)
      
      # if not installed from CRAN, install from Bioconductor
      if (!require(library, character.only = TRUE)) {
        suppressWarnings(source('http://bioconductor.org/biocLite.R')(library))
      }
    }
  }
}

#### set_global_variables ####
# Description: A function that sets all global variables for the aml project analysis. 
#
# Inputs: 
#     - locale = string representing where the analysis is being run
#
# Outputs: None
#
# Side-effects: 
#       Sets up global variables for the aml project analysis including the following: 
#         - RAW_DATA_DIRECTORY = folder where raw .fcs files are located
#         - DATA_DIRECTORY = folder where processed data files are located
#         - READER_OUTPUT = folder where figures from the reader functions should be saved
#         - HISTOGRAM_OUTPUT = folder where univariate histograms should be saved
set_global_variables <- function(locale = NULL){ 
  if (locale == "galaxia") { 
    
    CODE_DIRECTORY <<- here::here("scripts")
    RAW_DATA_DIRECTORY <<- file.path("~", "Box", "Tim", "Lab", "Data", "AML_Data", "AML_Compensated_Corrected_Myeloid")
    DATA_DIRECTORY <<- here::here("data")
    
    READER_OUTPUT <<- here::here("figures", "processing")
    
    HISTOGRAM_OUTPUT <<- here::here("figures", "histograms")
    
  } else if (locale == "sherlock") { 
    NULL
  } else { 
    NULL
  }
}

#marker_setup
#input: none
#output: none
#side-effects: Sets up global variables SURFACE_MARKERS, SIGNALING_MARKERS, TRX_FACTORS, and OTHER_MARKERS
#optimizations: set this up so that it can accept a string or an excel file with this information included

marker_setup <- function(){
  
  SURFACE_MARKERS <<- 
    c(
      "CD3_CD19",
      "CD10",
      "CD11b",
      "CD11c",
      "CD13",
      "CD14",
      "CD16",
      "CD33",
      "CD34",
      "CD38",
      "CD41",
      "CD45",
      "CD45RA",
      "CD47",
      "CD49f",
      "CD56",
      "CD61",
      "CD64",
      "CD68",
      "CD71",
      "CD90",
      "CD93",
      "CD99", 
      "CD109",
      "CD117", 
      "CD123",
      "CD135", 
      "CCR2", 
      "TIM-3", 
      #"PD-L1", 
      "HLA-DR"
    )
  
  SIGNALING_MARKERS <<- c("pAkt", "pCreb", "pErk", "pS6", "pSTAT3", "pSTAT5")
  
  TRX_FACTORS <<- c("CEBPa", "GATA-1", "PU.1")
  
  OTHER_MARKERS <<- c("caspase-3", "MPO")
  
  ALL_MARKERS <<- 
    c(
      SURFACE_MARKERS, 
      SIGNALING_MARKERS, 
      TRX_FACTORS, 
      OTHER_MARKERS
    )
  
  CLASSIFIER.MARKERS <<- 
    c(
      'CD45', 
      'CD34', 
      'CD38', 
      'CD61', 
      'CD14',
      'CD135', 
      'CD45RA', 
      'CD90', 
      'HLA-DR', 
      'CD41', 
      'CD13', 
      'CD11b', 
      'CD11c'
    )
  
  CLASSIFIER.POPULATIONS <<- 
    c(
      'HSC', 
      'MPP', 
      'CMP', 
      'GMP', 
      'MEP', 
      'Monocyte', 
      'DC', 
      'Macrophage',
      'Thrombocyte'
    )
}

#patient_setup.R
#input: none
#output: none
#side-effects: Sets up global variables HEALTHY CONTROLS, PAIRED_PATIENTS, UNPAIRED_PATIENTS
#optimizations: set this up so that it can accept a string or an excel file with this information included

patient_setup <- function(){
  
  HEALTHY_CONTROLS <<- 
    c(
      "bm5721", 
      "bm5871", 
      "bm6152", 
      "bm6631"
    )
  
  PAIRED_PATIENTS <<- 
    c(
      "parbiu", 
      "parwxu", 
      "parzuu", 
      "pasbee", 
      "pasgwh", 
      "patjsp", 
      "pasrtp", 
      "pasvzc", 
      "pasxvc", 
      "patgiy", 
      "pathiw", 
      "patlhb"
    )
  
  UNPAIRED_PATIENTS <<- 
    c(
      "papwhs", 
      "parajx", 
      "parant", 
      "parbfj", 
      "parcvp", 
      "parklc", 
      "parmme", 
      "parpwl", 
      "partxh", 
      "parxmp", 
      "parzwh", 
      "paseyc", 
      "paself", 
      "pastjj", 
      "paspsv", 
      "pasptm", 
      "pasxjy", 
      "pasyyw"
    )
}

#### source_tidyTOF.R ####
# Description: A function to source all prototype tidyTOF functions. 
#
# Inputs: 
#     - tidyTOF_directory = file path to the tidyTOF folder on your machine.
#
# Outputs: None
#
# Side-effects: All tidyTOF experimental functions are sourced into the currently active session.
source_tidyTOF <- function(tidyTOF_directory) {
  my_files <-
    list.files(tidyTOF_directory, pattern = "tof.+\\.R$", full.names = TRUE)
  
  for (file in my_files) { 
    source(file)
  }
}

