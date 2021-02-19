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
        suppressWarnings(BiocManager::install(library))
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
    
    HEALTHY_MYELOID_DIRECTORY <<- NULL
    
    HEALTHY_B_DIRECTORY <<- NULL
    
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
  
  CLASSIFIER_MARKERS <<- 
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
  
  CLASSIFIER_POPULATIONS <<- 
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
  
  ALL_CLASSIFIER_MARKERS <<- 
    c(
      "CD34", 
      "CD179a", 
      "CD38", 
      "IgMs", 
      "CD127", 
      "CD19", 
      "CD20", 
      "CD24", 
      "Tdt", 
      "IgMi", 
      "CD179b"
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

#### load_cancer_data() ####
# Description: A function to load data for the aml project
#
# Inputs: 
#     - cancer = string indicating if the data should be "aml" or "all" data
#     - raw = a boolean indicating if the data should be loaded from raw files or pre-loaded files 
#     - processing = string indicating the type of preprocessing to load the data with 
#     - sample = boolean indicating whether to load subsampled data or not
#
# Outputs: [m x n] tibble of m cells and n variables containing information about each cell
#
# Side-effects: 

load_cancer_data <- 
  function(
    cancer = NULL, 
    raw = FALSE, 
    processing = "arcsinh", 
    sample = FALSE
  ) { 
    
  # source tidyTOF functions 
  source_tidyTOF(file.path("~", "GitHub", "tidyTOF"))
  
  # handle processing options
  if (processing == "arcsinh") { 
    suffix <- "p" 
    my_transform <- function(x) asinh(x/5)
  } else if (all(processing == c("scale"))) { 
    suffix <- "s"
    my_transform <- scale
  } else if (all(processing == c("arcsinh", "scale"))) { 
    suffix <- "ps"
    my_transform <- 
      function(x) {
        asinh(x/5) %>% 
          scale() %>% 
          as.numeric()
      }
  } else if (all(processing == c("rank"))) { 
    suffix <- "_r"
    transform_fun <- rank
  }
  
  # handle raw or not raw (and resulting file paths)
  if (raw) { 
    if (str_to_lower(cancer) == "aml") {
      data_path <- 
        file.path(
          "~", 
          "Box", 
          "Tim", 
          "Lab", 
          "Data", 
          "AML_Data", 
          "AML_Compensated_Corrected_Myeloid"
        )
      
      result_tibble <-
        data_path %>% 
        tof_read_fcs() %>% 
        # rename and select out some variables that aren't needed
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
        rename(file_names = names) %>% 
        tof_preprocess(transform_fun = my_transform)
      
    } else if (str_to_lower(cancer) == "all") { 
      data_path <- 
        file.path(
          "~", 
          "Box", 
          "Tim", 
          "Lab", 
          "Data", 
          "DDPR_Data" 
        )
      
      marker_path <- here::here("docs", "ALL_panel.csv")
      
      marker_names <- 
        marker_path %>%
        read_csv() %>% 
        mutate(Metal = str_replace_all(Metal, "[()]", ""))
      
      my_data <- 
        data_path %>%
        list.files(path = ., full.names = TRUE) %>% 
        str_split(pattern = "_", simplify = TRUE) %>% 
        as_tibble() %>% 
        transmute(
          file_name = list.files(path = data_path, full.names = TRUE),
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
      
      result_tibble <- 
        my_data %>% 
        group_by(file_name) %>% 
        unnest(cols = data) %>% 
        ungroup() %>% 
        tof_preprocess(transform_fun = my_transform)
      
    } else { 
      stop("Only aml and all are valid options.")
    }
    
    # handle sampling 
    if (sample) {
      result_tibble <- 
        result_tibble %>% 
        group_by(file_names) %>% 
        slice_sample(n = 1000, replace = TRUE) %>% 
        ungroup()
    }
    
  } else {
    data_path <- 
      here::here("data") %>% 
      str_c("/", if_else(sample, "sampled_", ""), cancer, "_data_", suffix, ".rds")
    
    result_tibble <- 
      data_path %>% 
      read_rds()
  }
  
  # return result
  return(result_tibble)
  }


#### function_name() ####
# Description: 
#
# Inputs: 
#     - 
#     - 
#     -
#     - 
#
# Outputs:
#
# Side-effects: 
aml_iterate_ddpr_extraction <- function(cluster_col, threshold) { 
  
  result <- 
    tof_extract_ddpr(
      healthy_tibble = healthy_tibble, 
      cancer_tibble = cancer_tibble,
      cluster_col = {{cluster_col}}, 
      stimulation_col = stimulation, 
      lineage_cols = any_of(SURFACE_MARKERS), 
      signaling_cols = any_of(SIGNALING_MARKERS), 
      threshold = threshold
    )
  
  result
}

#### aml_extract_lineage_features() ####
# Description: 
#
# Inputs: 
#     - 
#     - 
#     -
#     - 
#
# Outputs:
#
# Side-effects: 
aml_extract_lineage_features <-
  function(cluster_col, central_tendency) {
    lineage_features <- 
      aml_data %>% 
      group_by(patient, condition) %>% 
      select(
        {{cluster_col}},
        any_of(c(SURFACE_MARKERS, TRX_FACTORS, OTHER_MARKERS))
      ) %>% 
      extract_feature_ct({{cluster_col}}, central_tendency = central_tendency) %>% 
      ungroup()
    
    basal_signal_features <- 
      aml_data %>% 
      filter(stimulation == "Basal") %>% 
      group_by(patient, condition) %>% 
      select(
        {{cluster_col}},
        any_of(SIGNALING_MARKERS)
      ) %>% 
      extract_feature_ct({{cluster_col}}, central_tendency = central_tendency) %>% 
      ungroup()
    
    result <- 
      lineage_features %>% 
      left_join(basal_signal_features, by = c("patient", "condition"))
    
    result
  }

#### function_name() ####
# Description: 
#
# Inputs: 
#     - 
#     - 
#     -
#     - 
#
# Outputs:
#
# Side-effects: 
aml_extract_signaling_features_emd <- 
  function(cluster_col) { 
    exports <- 
      c(
        "extract_feature_signaling", 
        "extract_feature_signaling_emd", 
        "tof_find_emd",
        "tof_find_density", 
        "SIGNALING_MARKERS"
      )
    
    emd_features <- 
      aml_data %>% 
      mutate(cluster_col = {{cluster_col}}) %>% 
      group_by(patient, condition) %>% 
      nest()
    
    # kind of hack-y, but will cause the extract_feature_signaling
    # function to look for the cluster_col column later instead of the 
    # original data variable, which foreach can't find (not sure why)
    cluster_col <- enquo(cluster_col)
    
    my_cluster <- parallel::makeCluster(14)
    registerDoParallel(my_cluster, 14)
    
    features <- 
      foreach(
        my_tibble = emd_features$data, 
        .packages = 
          c("dplyr", "purrr", "tidyr", "emdist", "rlang", "stringr"), 
        .export = exports
      ) %dopar% 
      extract_feature_signaling(
        tof_tibble = my_tibble,
        method = "emd",
        cluster_col = cluster_col,
        stimulation_col = stimulation,
        signaling_cols = all_of(SIGNALING_MARKERS),
        basal_level = "Basal"
      )
    
    parallel::stopCluster(my_cluster)
    
    emd_features <- 
      emd_features %>% 
      ungroup() %>% 
      transmute(patient, condition, features = features) %>% 
      unnest(cols = features)
    
    emd_features
  }


#### function_name() ####
# Description: 
#
# Inputs: 
#     - 
#     - 
#     -
#     - 
#
# Outputs:
#
# Side-effects: 
aml_extract_lineage_emd <- function(cluster_col, central_tendency) {
  
  result <- 
    aml_extract_lineage_features({{cluster_col}}, central_tendency) %>% 
    left_join(aml_extract_signaling_features_emd({{cluster_col}}))
  
  result
  
}

#### function_name() ####
# Description: 
#
# Inputs: 
#     - 
#     - 
#     -
#     - 
#
# Outputs:
#
# Side-effects: 
aml_extract_signaling_features_js <- function(cluster_col) { 
  
  result <- 
    aml_data %>% 
    group_by(patient, condition) %>% 
    extract_feature_signaling_js(
      cluster_col = {{cluster_col}}, 
      stimulation_col = stimulation, 
      signaling_cols = all_of(SIGNALING_MARKERS), 
      basal_level = "Basal"
    )
  
  result
}

#### function_name() ####
# Description: 
#
# Inputs: 
#     - 
#     - 
#     -
#     - 
#
# Outputs:
#
# Side-effects: 
aml_extract_lineage_js <- function(cluster_col, central_tendency) {
  
  result <- 
    aml_extract_lineage_features({{cluster_col}}, central_tendency) %>%
    left_join(aml_extract_signaling_features_js({{cluster_col}}))
  
  result
  
}

#### function_name() ####
# Description: 
#
# Inputs: 
#     - 
#     - 
#     -
#     - 
#
# Outputs:
#
# Side-effects:
aml_split_cv <- 
  function(feature_tibble, num_bootstraps) { 
    boots <- 
      feature_tibble %>% 
      bootstraps(times = num_bootstraps, strata = outcome)
    
    boots
  }


#### function_name() ####
# Description: 
#
# Inputs: 
#     - 
#     - 
#     -
#     - 
#
# Outputs:
#
# Side-effects:



#### function_name() ####
# Description: 
#
# Inputs: 
#     - 
#     - 
#     -
#     - 
#
# Outputs:
#
# Side-effects:
aml_create_recipe <- 
  function(feature_tibble, type = c("full", "baseline")) { 
    rlang::arg_match(type)
    
    if (type == "baseline") {
      aml_recipe <- 
        feature_tibble %>% 
        transmute(
          outcome = as.factor(outcome), 
          age, 
          wbc
        ) %>% 
        recipe(outcome ~ .) %>% 
        step_normalize(-outcome)
    }
    
    if (type == "full") {
      aml_recipe <- 
        feature_tibble %>% 
        mutate(outcome = as.factor(outcome)) %>% 
        recipe(formula = outcome ~ .) %>% 
        step_rm(
          patient, condition, 
          fab_category, cns_disease, 
          primary_cytogenetic_code, 
          cytogenetic_complexity
        ) %>% 
        step_normalize(-outcome)
    }
  }


#### function_name() ####
# Description: 
#
# Inputs: 
#     - 
#     - 
#     -
#     - 
#
# Outputs:
#
# Side-effects:

# CURRENTLY EXPERIMENTAL; DO NOT USE
aml_create_model <-
  function(feature_tibble, type = c("logistic_regression", "survival_regression")) { 
    rlang::arg_match(type)
    
    if (type == "logistic_regression") {
      ddpr_en_model_spec <- 
        logistic_reg(penalty = tune(), mixture = 1) %>%
        set_engine("glmnet") %>% 
        set_mode("classification")
      
    }
    
  }
