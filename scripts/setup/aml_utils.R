###     aml_utils.R
# Description: 
### Sources script files in the aml-cells project repository for a variety of 
### helper functions. 

# Author: Timothy Keyes
# Version: 2021-07-11

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



#flowSOM_cluster.R 

# Description: Performs FlowSOM clustering on an input dataset with a specified number k of nearest 
#              neighbors in the graph construction step. 
#
# Inputs:
#    expression_matrix = an [m x n] matrix of m cells and n measured signals on which the algorithm will operate. 
#    cluster_markers = character vector of the column (marker) names to use in the clustering.
#    num_cells_mc = number of cells from each cluster to use for monte carlo sampling and cluster metric calculations
#    make_patient_heatmap = boolean value indicating if a heatmap representing the percentage of each patient that falls in 
#                           each cluster should be plottted. Requires that expression_matrix has a column with the name 
#                           "patient." Default = FALSE.
#
# Returns: A flowCluster object with the following components... 
#     num_clusters = an integer indicating the number of clusters that were identified by the Phenograph algorithm
#     my_clusters: an integer vector of length m indicating the cluster assignment for each cell in expression_matrix
#     my_silhouettes: a double vector of length m indicating the silhouette value for each phenograph cluster
#     patient_heatmap: a ggplot object encoding a heatmap of the distribution of each patient's cells across 
#                      each of the phenograph clusters identified in the algorithmic step. If FALSE, will not return a heatmap. 
#     heatmap_matrix: the matrix that was used to construct the patient_heatmap.
#
#
# Optimizations: 
#
#     Currently, expression_matrix must contain a column called "patient" - change this to any categorical variable (or multiple).


flowSOM_cluster <- 
  function(
    expression_matrix = NULL, 
    cluster_markers = colnames(expression_matrix)
  ){
    
    my_FSO <- 
      list(
        data = data.matrix(expression_matrix), 
        compensate = FALSE, 
        spillover = NULL, 
        transform = FALSE, 
        scale = NULL, 
        prettyColnames = colnames(expression_matrix)
      )
    
    my_SOM <- 
      BuildSOM(
        fsom = my_FSO, 
        colsToUse = which(colnames(expression_matrix) %in% cluster_markers)
      )
    
    my_clusters <- my_SOM$map$mapping[,1] #flowSOM clustering must go here
    
    my_MST <- BuildMST(my_SOM, tSNE = FALSE)
    
    flowSOM_metaclusters <- 
      MetaClustering(
        data = my_MST$map$codes, 
        method = "metaClustering_consensus", 
        max = 50
      )
    
    flowSOM_metaclusters_cells <- flowSOM_metaclusters[my_MST$map$mapping[,1]]
    
    
    
    
    #deal with cells not assigned to a cluster - take out?
    if(length(my_clusters) != nrow(expression_matrix)) { 
      cell_string <- as.character(as.numeric(1:nrow(expression_matrix)))
      bad_cells <- which(!(cell_string %in% names(my_clusters)))
      bad_cell_clusters <- rep(NA, length(bad.cells)) 
      names(bad_cell_clusters) <- as.character(bad.cells)
      my_clusters <- c(my_clusters, bad_cell_clusters)
      my_clusters <- my_clusters[order(as.numeric(names(my_clusters)))] 
    }
    
    my_clusters <- 
      my_clusters %>% 
      as.numeric()
    
    #save num_clusters
    num_clusters <- length(unique(my_clusters))
    
    #return final list 
    flowCluster_result <- 
      list(
        num_clusters = num_clusters, 
        my_clusters = my_clusters, 
        my_metaclusters = flowSOM_metaclusters_cells
      )
    
    flowCluster_result
  } 

#
#
#
subset_sc_aml <- function(cluster_col, my_cluster, my_patient, my_condition, my_stimulation){ 
  aml_data %>% 
    filter(
      {{cluster_col}} == my_cluster,
      patient == my_patient,
      condition == my_condition, 
      stimulation %in% c("Basal", my_stimulation)
    )
}

#

# Description: 
#
# Inputs:
#    
#
# Returns: 
#
#
# Optimizations: 
#
explore_signaling_residuals <- 
  function(my_cluster_type, my_feature, my_outcome) {
  
  # obtain all data corresponding to the cluster_type you're interested in
  modeling_data <- 
    cancer_modeling_omnibus %>% 
    filter(cluster_type == my_cluster_type) %>% 
    pull(features) %>% 
    pluck(1) %>% 
    group_by(feature) %>% 
    nest()
  
  cluster_col <- 
    cancer_modeling_omnibus %>% 
    filter(cluster_type == my_cluster_type) %>% 
    pull(cluster_cols) %>% 
    pluck(1)
  
  # filter out only data corresponding to the feature and outcome of interest
  my_data <- 
    modeling_data %>% 
    filter(feature == my_feature) %>% 
    pull(data) %>% 
    pluck(1) %>% 
    select(patient, condition, ddpr, any_of(my_outcome)) %>% 
    drop_na()
  
  # model the ddpr features as a function of the my_outcome features
  linear_model <- 
    linear_reg() %>% 
    set_engine("lm")
  
  linear_recipe <- 
    recipe(formula = ddpr ~ ., x = my_data) %>% 
    step_rm(patient, condition) %>% 
    #step_scale(everything(), skip = TRUE) %>% 
    step_mutate_at(
      any_of(c("ddpr", my_outcome)), 
      fn = ~ .x / max(.x, na.rm = TRUE), 
      skip = TRUE
    ) %>% 
    step_naomit(everything(), skip = TRUE)
  
  linear_workflow <- 
    workflow() %>% 
    add_recipe(linear_recipe) %>% 
    add_model(linear_model)
  
  linear_fit <- 
    linear_workflow %>% 
    fit(data = my_data)
  
  my_residuals <- 
    linear_fit %>% 
    extract_model() %>% 
    pluck("residuals")
  
  my_predictions <- 
    linear_fit %>% 
    predict(new_data = my_data) %>% 
    pull(.pred)
  
  c(intercept, slope) %<-% 
    (
      linear_fit %>% 
        extract_model() %>% 
        pluck("coefficients") %>% 
        as.numeric() %>% 
        as.list()
    )
  
  modeling_result <- 
    linear_recipe %>%
    prep() %>% 
    juice() %>% 
    mutate(
      residuals = my_residuals, 
      predictions = my_predictions
    ) %>% 
    bind_cols(select(my_data, patient, condition))
  
  # make residual histogram
  residual_histogram <- 
    modeling_result %>% 
    ggplot(aes(x = residuals)) + 
    geom_histogram() + 
    labs(subtitle = my_feature)
  
  # make regression plot 
  regression_plot <- 
    modeling_result %>% 
    ggplot(aes_string(y = "ddpr", x = my_outcome)) + 
    geom_point() + 
    geom_abline(intercept = intercept, slope = slope, linetype = "dashed")
  
  # find all the variables we need to filter the single-cell data
  c(my_stim, my_marker, my_cluster) %<-% 
    str_split(my_feature, pattern = "_", simplify = TRUE)
  
  aml_subset <- 
    modeling_result %>% 
    mutate(
      data = 
        map2(
          .x = patient,
          .y = condition,
          .f = ~ 
            subset_sc_aml(
              cluster_col = {{cluster_col}}, 
              my_cluster = my_cluster, 
              my_patient = .x, 
              my_condition = .y, 
              my_stimulation = my_stim
            ) %>% 
            select(patient, condition, stimulation, {{cluster_col}}, any_of(my_marker))
        ), 
      num_cells = map_int(
        .x = data, 
        .f = ~ 
          count(.x, stimulation) %>% 
          pull(n) %>% 
          min()
      )
    )
  
  correlation <- 
    aml_subset %>% 
    summarize(correlation = cor(residuals, num_cells)) %>% 
    pull(correlation) %>% 
    round(3)
  
  # make a scatterplot of number of cells vs. the residuals
  correlation_scatterplot <- 
    aml_subset %>%
    ggplot(aes(x = num_cells, y = residuals)) + 
    geom_point() + 
    labs(caption = str_c("Correlation = ", correlation))
    
  correlation_scatterplot <- 
    correlation_scatterplot %>% 
    ggMarginal(type = "histogram")
  
  # make histograms comparing the highest, lowest,and median residual distributions 
  # with one another
  my_length <- length(aml_subset$residuals)
  
  max_residual <- max(abs(aml_subset$residuals))
  min_residual <- min(abs(aml_subset$residuals))
  mid_residual <- 
    abs(aml_subset$residuals[rank(abs(aml_subset$residuals)) == ceiling(my_length / 2)])
  
  patient_order <- 
    aml_subset %>% 
    filter(abs(residuals) %in% c(max_residual, min_residual, mid_residual)) %>% 
    arrange(-abs(residuals)) %>% 
    pull(patient)
  
  annotation_tibble <- 
    aml_subset %>% 
    filter(abs(residuals) %in% c(max_residual, min_residual, mid_residual)) %>% 
    select(data) %>% 
    unnest(cols = data) %>% 
    mutate(patient = factor(patient, levels = patient_order)) %>% 
    group_by(patient, condition, stimulation) %>% 
    summarize(
      num_cells = n(), 
      across(any_of(my_marker), ~ quantile(.x, 0.85))
    ) %>% 
    rename_with(.fn = ~ return("channel"), .cols = any_of(my_marker))


  density_plot <- 
    aml_subset %>%
    filter(abs(residuals) %in% c(max_residual, min_residual, mid_residual)) %>% 
    mutate(patient = factor(patient, levels = patient_order)) %>% 
    unnest() %>%
    ggplot(aes_string(x = my_marker, fill = "stimulation")) + 
    geom_vline(xintercept = asinh(10/5), linetype = "dashed", color = "black") + 
    geom_density(aes(y = ..scaled..)) + 
    geom_text(
      aes(x = asinh(10/5), y = 1, label = str_c(num_cells, " cells")), 
      data = annotation_tibble, 
      hjust = 0, 
      nudge_x = 0.1
    ) + 
    facet_grid(rows = vars(stimulation), cols = vars(patient), scales = "free") + 
    labs(
      subtitle = str_c(my_marker, " in this cluster: ", my_cluster),
      y = "Scaled density",
      caption = "From right to left: Maximum residual, median residual, minimum residual"
    )
  
  final_result <- 
    tibble(
      cluster_type = my_cluster_type, 
      feature = my_feature,
      modeling_result = list(modeling_result),
      regression_plot = list(regression_plot),
      residual_histogram = list(residual_histogram), 
      correlation_scatterplot = list(correlation_scatterplot), 
      density_plot = list(density_plot)
    )
  return(final_result)
}


