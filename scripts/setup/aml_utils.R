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

#### Another function ####

