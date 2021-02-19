###     feature_extraction_utils.R
# Description: 
### Creates some functions for feature extraction 

# Author: Timothy Keyes
# Version: 2020-11-05



# Description: A function to extract surface-level features as described in the DDPR paper.  
#
# Inputs: 
#     -tof_tibble = a single-cell data frame 
#     -lineage_markers = tidy selection of column names of lineage markers (surface markers and
#                        transcription factors)
#     -central_tendency = a function indicating the measure of central tendency you want to measure 
#                         in each cluster for each of your features 
#     -patient_var = unquoted name of the column encoding each patient's id
#     -cluster_var = unquoted name of the column encoding the cluster id to which each cell belongs
#     -expanded_clusters = character vector of cluster ids indicating which clusters are expanded 
#                          relative to healthy samples. 
#
# Outputs: lineage_tibble = a tibble containing the patient-level lineage features for the input dataset.
#
# Side-effects: None
#
# Notes: This function is designed to be easy to use for anyone who isn't used to tidy syntax, but it 
#        also adheres to tidy rule. For example, you can group before using this function to extract features 
#        by additional grouping variables (in addition to patient and cluster). 
extract_ddpr_lineage <- 
  function(
    tof_tibble, # entire single-cell data frame
    lineage_markers, # tidy selection of column names of lineage markers (surface and trx factors)
    central_tendency = mean,
    patient_var = patient, 
    cluster_var = cluster, 
    expanded_clusters # character vector of clusters that were expanded from the abundance calculation
  ) { 
    lineage_tibble <- 
      tof_tibble %>% 
      filter({{cluster_var}} %in% expanded_clusters) %>% 
      group_by({{patient_var}}, {{cluster_var}}, .add = TRUE) %>% 
      summarize(
        across({{lineage_markers}}, central_tendency), 
        .groups = "drop_last"
      ) %>% 
      pivot_longer(
        cols = {{lineage_markers}}, 
        names_to = "marker", 
        values_to = "expression"
      ) %>% 
      pivot_wider(
        names_from = c({{cluster_var}}, marker), 
        values_from = expression, 
        names_sep = "_"
      ) %>% 
      ungroup()
    return(lineage_tibble)
  }


# Description: A function to extract signaling features as described in the DDPR paper. 
#
# Inputs: 
#     -tof_tibble = a single-cell data frame 
#     -signaling_markers = tidy selection of column names of signaling markers
#     -central_tendency = a function indicating the measure of central tendency you want to measure 
#                         in each cluster for each of your features 
#     -patient_var = unquoted name of the column encoding each patient's id
#     -cluster_var = unquoted name of the column encoding the cluster id to which each cell belongs
#     -stimulation_var = unquoted name of the column encoding the stimulation used for each cell
#     -threshold = constant threshold above which a cell's signaling should be considered "positive"
#     -stimulation_baseline = character vector of values within the stimulation_column that denote what the "basal" condition
#                             is called in your dataset. 
#     -expanded_clusters = character vector of cluster ids indicating which clusters are expanded 
#                          relative to healthy samples. 
#
# Outputs: signaling_tibble = a tibble containing the patient-level signaling features for the input dataset.
#
# Side-effects: None
#
# Notes: This function is designed to be easy to use for anyone who isn't used to tidy syntax, but it 
#        also adheres to tidy rule. For example, you can group before using this function to extract features 
#        by additional grouping variables (in addition to patient and cluster). 

extract_ddpr_signaling <- 
  function(
    tof_tibble, 
    signaling_markers, 
    patient_var = patient, 
    cluster_var = cluster,
    stimulation_var = stimulation,
    threshold = asinh(10 / 5),
    stimulation_baseline = c("Basal", "basal"), 
    expanded_clusters 
  ) {
    signaling_tibble <- 
      tof_tibble %>% 
      filter({{cluster_var}} %in% expanded_clusters) %>% 
      group_by({{patient_var}}, {{cluster_var}}, {{stimulation_var}}, .add = TRUE) %>% 
      summarize(
        across(all_of(signaling_markers), ~ mean(.x > threshold)), 
        .groups = "drop_last"
      ) %>% 
      pivot_longer(
        cols = all_of(signaling_markers), 
        names_to = "marker",
        values_to = "prop_positive"
      ) %>% 
      pivot_wider(
        names_from = c({{cluster_var}}, {{stimulation_var}}, marker), 
        values_from = prop_positive, 
        names_sep = "_"
      ) %>% 
      ungroup()
    
    return(signaling_tibble)
    
  }


# Description: A function to extract signaling features as described in the DDPR paper. 
#
# Inputs: 
#     -tof_tibble = a single-cell data frame 
#     -patient_var = unquoted name of the column encoding each patient's id
#     -cluster_var = unquoted name of the column encoding the cluster id to which each cell belongs
#
# Outputs: abundance_tibble = a tibble containing the proportion of each cell cluster in each patient.
#
# Side-effects: None
#
# Notes: This function is designed to be easy to use for anyone who isn't used to tidy syntax, but it 
#        also adheres to tidy rule. For example, you can group before using this function to extract features 
#        by additional grouping variables (in addition to patient and cluster). 
extract_cluster_abundances <- 
  function(
    tof_tibble,
    patient_var = patient, 
    cluster_var = cluster
  ) { 
    abundance_tibble <- 
      tof_tibble %>% 
      count({{patient_var}}, {{cluster_var}}, name = "abundance") %>% 
      group_by({{patient_var}}) %>% 
      mutate(
        abundance = abundance / sum(abundance), 
        "{{cluster_var}}" := str_c({{cluster_var}}, "abundance", sep = "_")
      ) %>% 
      pivot_wider(
        names_from = {{cluster_var}},  
        values_from = abundance
      ) %>% 
      ungroup()
    
    return(abundance_tibble)
    
  }

# Description: 
### Performs feature extraction from a tof_tibble in prep for predictive modeling according to the DDPR procedure.
#
# Inputs: 
#     - tof_tibble = a tibble, data.frame, or something that can be coerced into either 
#       of these
#     - method = what feature matrix do you want to extract? Currently implemented: "ddpr", "citrus", 
#                "js", and "emd" 
#     - lineage_markers = A vector of non-quoted variables representing columns that contain
#                      single-cell protein measurements corresponding to phenotypic/lineage markers
#                      that *do not* differ between stimulation conditions. Anything that works in the 
#                      first argument of dplyr::across will work. See ?across. 
#                      Supports tidy selection using tidy_select helpers. The default is to 
#                      transform all numeric columns. 
#     - signaling_markers = A vector of non-quoted variables representing columns that contain
#                      single-cell protein measurements corresponding to functional/signaling markers 
#                      that *will* differ between stimulation conditions. Anything that works in the 
#                      first argument of dplyr::across will work. See ?across. 
#                      Supports tidy selection using tidy_select helpers. The default is to 
#                      transform all numeric columns. 
#     - patient_var = unquoted name of the variable encoding the patient IDs.
#     - condition_var = unquoted name of the variable encoding the condition (healthy or not) labels for each cell.                         each cell.
#     - cluster_var = unquoted name of the variable encoding the cluster labels.
#     - stimulation_var = unquoted name of the variable encoding the stimulation labels for each cell.
#     - condition_healthy = a string representing the healthy value in the stimulation column. 
#     - stimulation_baseline = a string representing the baseline value in the stimulation column. 
#
# Outputs: 
#     - feature_tibble = a tibble containing the calculated features for each patient
#
# Dependencies: 
#     - tidyverse library
#
#
tof_feature_extract_ddpr <- 
  function(
    tof_tibble, 
    lineage_markers, 
    signaling_markers, 
    patient_var = patient,
    condition_var = condition,
    cluster_var = NULL,
    stimulation_var = stimulation, 
    condition_healthy = c("healthy", "Healthy"),
    stimulation_baseline = c("basal", "Basal"),
    central_tendency = mean,
    threshold = asinh(10 / 5),
    ...
  ){ 
    abundance_features <- 
      tof_tibble %>% 
      group_by({{condition_var}}) %>% 
      extract_cluster_abundances(
        patient_var = {{patient_var}}, 
        cluster_var = {{cluster_var}}
      ) %>%
      mutate(across(everything(), replace_na, replace = 0))
    
    #find expanded clusters
    expanded_clusters <- 
      abundance_features %>% 
      group_by(
        is_healthy = 
          if_else(
            {{condition_var}} %in% condition_healthy, 
            "healthy", 
            "not_healthy"
          )
      ) %>% 
      summarize(across(c(-{{condition_var}}, -{{patient_var}}), mean)) %>% 
      pivot_longer(
        cols = -is_healthy, 
        names_to = "cluster", 
        values_to = "proportion" 
      ) %>% 
      pivot_wider(
        names_from = is_healthy,
        values_from = proportion
      ) %>% 
      filter(healthy < not_healthy) %>% 
      mutate(cluster = str_remove(cluster, "_abundance")) %>% 
      pull(cluster)
    
    # find lineage features
    lineage_features <- 
      tof_tibble %>% 
      extract_ddpr_lineage(
        lineage_markers = {{lineage_markers}}, #may have to unquote this
        central_tendency = central_tendency,
        patient_var = {{patient_var}}, 
        cluster_var = {{cluster_var}}, 
        expanded_clusters = expanded_clusters
      )
    
    # find signaling features
    signaling_features <- 
      tof_tibble %>% 
      extract_ddpr_signaling(
        signaling_markers = {{signaling_markers}}, 
        patient_var = {{patient_var}}, 
        cluster_var = {{cluster_var}},
        stimulation_var = {{stimulation_var}},
        threshold = threshold,
        stimulation_baseline = stimulation_baseline, 
        expanded_clusters = expanded_clusters
      )
    
    # combine all features
    feature_tibble <- 
      abundance_features %>% 
      left_join(lineage_features) %>% 
      left_join(signaling_features)
    
    return(feature_tibble)
  }



# Description: 
### Performs feature extraction from a tof_tibble in prep for predictive modeling according to the CITRUS procedure.
#
# Inputs: 
#     - tof_tibble = a tibble, data.frame, or something that can be coerced into either 
#       of these
#     - method = what feature matrix do you want to extract? Currently implemented: "ddpr", "citrus", 
#                "js", and "emd" 
#     - lineage_markers = A vector of non-quoted variables representing columns that contain
#                      single-cell protein measurements corresponding to phenotypic/lineage markers
#                      that *do not* differ between stimulation conditions. Anything that works in the 
#                      first argument of dplyr::across will work. See ?across. 
#                      Supports tidy selection using tidy_select helpers. The default is to 
#                      transform all numeric columns. 
#     - signaling_markers = A vector of non-quoted variables representing columns that contain
#                      single-cell protein measurements corresponding to functional/signaling markers 
#                      that *will* differ between stimulation conditions. Anything that works in the 
#                      first argument of dplyr::across will work. See ?across. 
#                      Supports tidy selection using tidy_select helpers. The default is to 
#                      transform all numeric columns. 
#     - patient_var = unquoted name of the variable encoding the patient IDs.
#     - condition_var = unquoted name of the variable encoding the condition (healthy or not) labels for each cell.                         each cell.
#     - cluster_var = unquoted name of the variable encoding the cluster labels.
#     - stimulation_var = unquoted name of the variable encoding the stimulation labels for each cell.
#     - condition_healthy = a string representing the healthy value in the stimulation column. 
#     - stimulation_baseline = a string representing the baseline value in the stimulation column. 
#
# Outputs: 
#     - feature_tibble = a tibble containing the calculated features for each patient
#
# Dependencies: 
#     - tidyverse library
#
#
tof_feature_extract_citrus <- 
  function(
    tof_tibble, 
    lineage_markers, 
    signaling_markers, 
    patient_var = patient,
    condition_var = condition,
    cluster_var = NULL,
    stimulation_var = stimulation, 
    condition_healthy = c("healthy", "Healthy"),
    stimulation_baseline = c("basal", "Basal"),
    central_tendency = mean,
    threshold = asinh(10 / 5),
    ...
  ){ 
    
  }


# Description: 
### Performs feature extraction from a tof_tibble in prep for predictive modeling.   
#
# Inputs: 
#     - tof_tibble = a tibble, data.frame, or something that can be coerced into either 
#       of these
#     - method = what feature matrix do you want to extract? Currently implemented: "ddpr", "citrus", 
#                "js", and "emd" 
#     - lineage_markers = A vector of non-quoted variables representing columns that contain
#                      single-cell protein measurements corresponding to phenotypic/lineage markers
#                      that *do not* differ between stimulation conditions. Anything that works in the 
#                      first argument of dplyr::across will work. See ?across. 
#                      Supports tidy selection using tidy_select helpers. The default is to 
#                      transform all numeric columns. 
#     - signaling_markers = A vector of non-quoted variables representing columns that contain
#                      single-cell protein measurements corresponding to functional/signaling markers 
#                      that *will* differ between stimulation conditions. Anything that works in the 
#                      first argument of dplyr::across will work. See ?across. 
#                      Supports tidy selection using tidy_select helpers. The default is to 
#                      transform all numeric columns. 
#     - patient_var = unquoted name of the variable encoding the patient IDs.
#     - condition_var = unquoted name of the variable encoding the condition (healthy or not) labels for each cell.                         each cell.
#     - cluster_var = unquoted name of the variable encoding the cluster labels.
#     - stimulation_var = unquoted name of the variable encoding the stimulation labels for each cell.
#     - condition_healthy = a string representing the healthy value in the stimulation column. 
#     - stimulation_baseline = a string representing the baseline value in the stimulation column. 
#
# Outputs: 
#     - feature_tibble = a tibble containing the calculated features for each patient
#
# Dependencies: 
#     - tidyverse library
#
#
tof_feature_extract <- 
  function(
    tof_tibble, 
    method, 
    lineage_markers, 
    signaling_markers, 
    patient_var = patient,
    condition_var = condition,
    cluster_var = NULL,
    stimulation_var = stimulation, 
    condition_healthy = c("healthy", "Healthy"),
    stimulation_baseline = c("basal", "Basal"),
    central_tendency = mean,
    threshold = asinh(10 / 5),
    ...
  ){ 
    
    # setup stuff 
    
    
    
    # method-specific calculations
    
    if (method == "ddpr") { 
      feature_tibble <-
        tof_tibble %>% 
        tof_feature_extract_ddpr(
          lineage_markers = {{lineage_markers}}, 
          signaling_markers = {{signaling_markers}}, 
          patient_var = {{patient_var}},
          condition_var = {{condition_var}},
          cluster_var = {{cluster_var}},
          stimulation_var = {{stimulation_var}}, 
          condition_healthy = condition_healthy,
          stimulation_baseline = stimulation_baseline,
          central_tendency = central_tendency,
          threshold = threshold
        )
    
    } else if (method == "citrus") { 
      
      
    } else if (method == "js") { 
      NULL
    } else if (method == "") {
      NULL
    }
    
    return(feature_tibble)
    
  }


# for testing

temp <- 
  aml_data %>% 
  tof_feature_extract(
    method = "ddpr", 
    lineage_markers = all_of(SURFACE_MARKERS), 
    signaling_markers = all_of(SIGNALING_MARKERS), 
    cluster_var = mahalanobis_cluster 
  )

