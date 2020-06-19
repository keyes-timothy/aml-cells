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
#         - 
set_global_variables <- function(locale = NULL){ 
  if (locale == "galaxia") { 
    RAW_DATA_DIRECTORY <<- file.path("~", "Box", "Tim", "Lab", "Data", "AML_Data", "AML_Compensated_Corrected_Myeloid")
    DATA_DIRECTORY <<- here::here("data")
    
    READER_OUTPUT <<- here::here("figures", "processing")
    
  } else if (locale == "sherlock") { 
    NULL
  } else { 
    NULL
  }
}

  

