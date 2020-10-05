AML Feature Extraction Vignette
================
tkeyes
2020-10-05

  - [Executive Summary](#executive-summary)
  - [Read in data](#read-in-data)
  - [Analysis](#analysis)
      - [DDPR feature extraction](#ddpr-feature-extraction)
  - [Developmental Classifier
    matrices](#developmental-classifier-matrices)
      - [DDPR extraction](#ddpr-extraction)
      - [CITRUS Extraction](#citrus-extraction)
      - [Other Extraction](#other-extraction)

``` r
# Libraries
libraries <- 
  c(
    "flowCore", 
    "tidyverse", 
    "readxl", 
    "ggridges", 
    "rlang", 
    "ggthemes", 
    "foreach", 
    "doParallel", 
    "ggrepel", 
    "ggiraph", 
    "ggiraphExtra", 
    "FlowSOM"
  )

source(here::here("scripts", "setup", "aml_utils.R")) #may need to change on scg machines
call_libraries(libraries)

# Parameters
set_global_variables(locale = "galaxia")
md_path <- here::here("data-raw", "AML_metadata.xlsx")
tidyTOF_directory <- file.path("~", "GitHub", "tidyTOF")

CLUSTER_OUTPUT <- here::here("data", "cancer_data_clustered.rds")

# Sourcing tidyTOF functions
source_tidyTOF(tidyTOF_directory)

# set up aml marker and patient information
marker_setup()
patient_setup()

# Misc globals for running different parts of this vignette
is_sampled <- FALSE
aml_path <- here::here("data", "aml_data_clustered.rds")
all_path <- here::here("data", "all_data_clustered.rds")

if (is_sampled) { 
  aml_path <- aml_path %>% str_replace(".rds", "_sampled.rds")
  all_path <- all_path %>% str_replace(".rds", "_sampled.rds")
}
```

# Executive Summary

In this document, we formalize the feature extraction methods used in a
variety of CyTOF papers as a preamble to using them in some predictive
models. The feature extraction methods that we use include the
following:

  - [Developmentally-Dependent Predictor of Relapse
    (DDPR)](https://pubmed.ncbi.nlm.nih.gov/29505032/)
  - [CITRUS](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4084463/)/[Statistical
    Scaffold](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4537647/)
  - Some hybrid methods using [Earth-Mover’s
    Distance](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4805242/) or
    the [Jensen-Shannon
    divergence](https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence)
    for signaling features (which was previously used in [this
    paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4076922/))

# Read in data

``` r
aml_data <- 
  aml_path %>% 
  read_rds()

all_data <- 
  all_path %>% 
  read_rds()
```

# Analysis

``` r
aml_data %>% colnames()
```

    ##  [1] "CD45"                        "CD61"                       
    ##  [3] "CD99"                        "127I"                       
    ##  [5] "CD45RA"                      "CD93"                       
    ##  [7] "CD3_CD19"                    "CCR2"                       
    ##  [9] "CD117"                       "CD123"                      
    ## [11] "CD64"                        "CD90"                       
    ## [13] "CD38"                        "CD34"                       
    ## [15] "CEBPa"                       "pSTAT5"                     
    ## [17] "CD11c"                       "CD13"                       
    ## [19] "pAkt"                        "TIM-3"                      
    ## [21] "CD56"                        "CD10"                       
    ## [23] "PU.1"                        "CD33"                       
    ## [25] "CD14"                        "caspase-3"                  
    ## [27] "GATA-1"                      "pSTAT3"                     
    ## [29] "CD41"                        "CD16"                       
    ## [31] "CD68"                        "MPO"                        
    ## [33] "pErk"                        "CD47"                       
    ## [35] "CD135"                       "CD109"                      
    ## [37] "pS6"                         "CD49f"                      
    ## [39] "HLA-DR"                      "CD71"                       
    ## [41] "pCreb"                       "191Ir"                      
    ## [43] "193Ir"                       "cisplatin"                  
    ## [45] "CD11b"                       "file_names"                 
    ## [47] "plate"                       "patient"                    
    ## [49] "stimulation"                 "condition"                  
    ## [51] "cell_id"                     "cluster_flowSOM_surface"    
    ## [53] "cluster_flowSOM_signaling"   "cluster_flowSOM_mah"        
    ## [55] "cluster_flowSOM_cos"         "cluster_flowSOM_all_markers"
    ## [57] "mahalanobis_HSC"             "mahalanobis_MPP"            
    ## [59] "mahalanobis_CMP"             "mahalanobis_GMP"            
    ## [61] "mahalanobis_MEP"             "mahalanobis_Monocyte"       
    ## [63] "mahalanobis_DC"              "mahalanobis_Macrophage"     
    ## [65] "mahalanobis_Thrombocyte"     "mahalanobis_cluster"        
    ## [67] "cosine_HSC"                  "cosine_MPP"                 
    ## [69] "cosine_CMP"                  "cosine_GMP"                 
    ## [71] "cosine_MEP"                  "cosine_Monocyte"            
    ## [73] "cosine_DC"                   "cosine_Macrophage"          
    ## [75] "cosine_Thrombocyte"          "cosine_cluster"

``` r
# DDPR helper functions
central_tendency <- mean
ddpr_surface_extract <- 
  function(
    patient_tibble, 
    central_tendency = mean, 
    cluster_var = cluster
  ) { 
    surface_tibble <- 
      patient_tibble %>% 
      group_by({{cluster_var}}) %>% 
       summarize(across(everything(), central_tendency), .groups = "drop_last") %>% 
       pivot_longer(
         cols = -{{cluster_var}}, 
         names_to = "marker", 
         values_to = "expression"
       ) %>% 
       pivot_wider(
         names_from = c({{cluster_var}}, marker), 
         values_from = expression, 
         names_sep = "_"
       )
    return(surface_tibble)
  }


# test to make sure ... works in tof_extract_feature_matrix()
ddpr_signaling_extract <- 
  function(
    patient_tibble, 
    central_tendency = mean, 
    cluster_var = cluster,
    stimulation_var = stimulation,
    threshold = asinh(10 / 5), 
    stimulation_baseline = c("Basal", "basal")
  ) {
    signaling_tibble <- 
      patient_tibble %>% 
      filter(!({{stimulation_var}} %in% stimulation_baseline)) %>% 
      mutate(
        cluster_stimulation = 
          str_c({{cluster_var}}, {{stimulation_var}}, sep = "_")
      ) %>% 
      select(-{{stimulation_var}}, -{{cluster_var}}) %>% 
      group_by(cluster_stimulation) %>% 
      summarize(
        across(everything(), ~ mean(.x > threshold)), 
        .groups = "drop_last"
      ) %>% 
      pivot_longer(
        cols = -cluster_stimulation, 
        names_to = "marker",
        values_to = "prop_positive"
      ) %>% 
      pivot_wider(
        names_from = c(cluster_stimulation, marker), 
        values_from = prop_positive, 
        names_sep = "_"
      ) 
    
    basal_values <- 
      patient_tibble %>% 
      filter({{stimulation_var}} %in% stimulation_baseline) %>% 
      mutate(
        cluster_stimulation = 
          str_c({{cluster_var}}, {{stimulation_var}}, sep = "_")
      ) %>% 
      select(-{{stimulation_var}}, -{{cluster_var}}) %>% 
      group_by(cluster_stimulation) %>% 
      summarize(
        across(everything(), central_tendency), 
        .groups = "drop_last"
      ) %>% 
      pivot_longer(
        cols = -cluster_stimulation, 
        names_to = "marker", 
        values_to = "expression"
      ) %>% 
      pivot_wider(
        names_from = c(cluster_stimulation, marker), 
        values_from = expression, 
        names_sep = "_"
      )
    
    signaling_tibble <- bind_cols(basal_values, signaling_tibble)
    
    return(signaling_tibble)
      
  }


# CITRUS helper functions - needs to be adjusted to find the difference between basal and each stim. 

citrus_feature_extract <- 
  function(
    patient_tibble, 
    central_tendency = mean, 
    cluster_var = cluster
  ) {
    citrus_tibble <- 
      patient_tibble %>% 
      group_by({{cluster_var}}) %>% 
      summarize(across(everything(), central_tendency), .groups = "drop_last") %>% 
      pivot_longer(
        cols = -{{cluster_var}}, 
        names_to = "marker", 
        values_to = "expression"
      ) %>% 
      pivot_wider(
        names_from = c({{cluster_var}}, marker), 
        values_from = expression, 
        names_sep = "_"
      )
    return(citrus_tibble)
    
  }




extract_cluster_abundances <- 
  function(
    patient_tibble, 
    cluster_var = cluster
  ) { 
    abundance_vector <- 
      patient_tibble %>% 
      count({{cluster_var}}, name = "abundance") %>% 
      mutate(
        abundance = abundance / sum(abundance), 
        "{{cluster_var}}" := str_c({{cluster_var}}, "abundance", sep = "_")
      ) %>% 
      pivot_wider(
        names_from = {{cluster_var}},  
        values_from = abundance
      )
    
    return(abundance_vector)
    
  }
```

``` r
# DDPR - have to replace mahalanobis cluster with general cluster_var when making function

feature_tibble <- 
  aml_data %>% 
  group_by(patient) %>% 
  slice_sample(prop = 0.1) %>% 
  ungroup() %>% 
  select(
    all_of(c(SURFACE_MARKERS, TRX_FACTORS, SIGNALING_MARKERS, "MPO")),
    patient, 
    stimulation, 
    condition, 
    mahalanobis_cluster
  ) %>% 
  nest(
    surface_data = c(all_of(c(SURFACE_MARKERS, TRX_FACTORS, "MPO")), mahalanobis_cluster), 
    signaling_data = c(all_of(SIGNALING_MARKERS), mahalanobis_cluster, stimulation),
  ) %>% 
  transmute(
    patient, 
    condition, 
    surface_vector = 
      map(
        .x = surface_data, 
        .f = ddpr_surface_extract, 
        cluster_var = mahalanobis_cluster
      ), 
    signaling_vector = 
      map(
        .x = signaling_data, 
        .f = ddpr_signaling_extract, 
        cluster_var = mahalanobis_cluster
      ), 
    abundance_vector = 
      map(
        .x = surface_data, 
        .f = extract_cluster_abundances, 
        cluster_var = mahalanobis_cluster
      )
  ) %>% 
  unnest(surface_vector, signaling_vector, abundance_vector)
```

## DDPR feature extraction

In this approach, the feature matrix is extracted from the single-cell
dataset with a few rules:

  - For surface markers and transcription factors, calculate the \_\_\_
    for each protein within each cluster within each patient. Do this
    across all cells collected from each patient regardless of what
    stimulation condition they were subject to.
  - For phospho-signals, calculate the percentage of “positive” cells in
    each population for each stimulation condition using a simple
    threshold (in the original paper, the threshold was 10 counts).
  - In addition to what’s described above, for the “basal” stimulation
    condition, calculate the \_\_\_ for each phospho-signal and include
    that as a feature as well.

<!-- end list -->

``` r
# Description: 
### Performs feature extraction from a tof_tibble in prep for predictive modeling.   
#
# Inputs: 
#     - tof_tibble = a tibble, data.frame, or somehing that can be coerced into either 
#       of these
#     - method = what feature matrix do you want to extract? Currently implemented: "ddpr", "citrus", 
#                "js", and "emd" 
#     - phenotypic_markers = A vector of non-quoted variables representing columns that contain
#                      single-cell protein measurements corresponding to phenotypic/lineage markers
#                      that *do not* differ between stimnulation conditions. Anything that works in the 
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
#     - cluster_var = unquoted name of the variable encoding the cluster labels.
#     - stimulation_var = unquoted name of the variable encoding the stimulation labels for each cell.
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
    phenotypic_markers, 
    signaling_markers, 
    patient_var = patient,
    condition_var = condition,
    cluster_var = NULL,
    stimulation_var = stimulation, 
    stimulation_baseline = c("basal", "Basal"),
    ...
  ){ 
    
    # setup stuff 
    
    
    
    # method-specific calculations
    
    if (method == "ddpr") { 
      feature_tibble <- 
        tof_tibble %>% 
        select(
          {{phenotypic_markers}},
          {{signaling_markers}}, 
          {{patient_var}}, 
          {{stimulation_var}}, 
          {{condition_var}},
          {{cluster_var}}
        ) %>% 
        nest(
          surface_data = c({{phenotypic_markers}}, {{cluster_var}}), 
          signaling_data = c({{signaling_markers}}, {{cluster_var}}, {{stimulation_var}}),
        ) %>% 
        transmute(
          {{patient_var}}, 
          {{condition_var}}, 
          surface_vector = 
            map(
              .x = surface_data, 
              .f = ddpr_surface_extract, 
              cluster_var = {{cluster_var}}
            ), 
          signaling_vector = 
            map(
              .x = signaling_data, 
              .f = ddpr_signaling_extract, 
              cluster_var = {{cluster_var}}
            ), 
          abundance_vector = 
            map(
              .x = surface_data, 
              .f = extract_cluster_abundances, 
              cluster_var = {{cluster_var}}
            )
        ) %>% 
        unnest(c(surface_vector, signaling_vector, abundance_vector))
      
      return(feature_tibble)
      
    } else if (method == "citrus") { 
      
      feature_tibble <- 
        tof_tibble %>% 
        select(
          {{phenotypic_markers}},
          {{signaling_markers}}, 
          {{patient_var}}, 
          {{stimulation_var}}, 
          {{condition_var}},
          {{cluster_var}}
        ) %>% 
        nest(
          marker_data = 
            c(
              {{phenotypic_markers}}, 
              {{signaling_markers}}, 
              {{cluster_var}}
            )
        ) %>% 
        transmute(
          {{patient_var}}, 
          {{condition_var}}, 
          data_vector = 
            map(
              .x = marker_data, 
              .f = citrus_feature_extract, 
              cluster_var = {{cluster_var}}
            ), 
          abundance_vector = 
            map(
              .x = marker_data, 
              .f = extract_cluster_abundances, 
              cluster_var = {{cluster_var}}
            )
        ) %>% 
        unnest(c(data_vector, abundance_vector))
      
    } else if (method == "js") { 
      NULL
    } else if (method == "") {
      NULL
    }
    
  }
```

# Developmental Classifier matrices

## DDPR extraction

``` r
mahalanobis_dev_matrix_10 <- 
  aml_data %>% 
  tof_feature_extract(
    method = "ddpr", 
    phenotypic_markers = c(all_of(c(SURFACE_MARKERS, TRX_FACTORS, "MPO"))), 
    signaling_markers = all_of(SIGNALING_MARKERS), 
    patient_var = patient, 
    condition_var = condition,
    cluster_var = mahalanobis_cluster, 
    stimulation_var = stimulation, 
    stimulation_baseline = "Basal"
  )

write_rds(
  x = mahalanobis_dev_matrix_10, 
  path = here::here("data", "mahalanobis_dev_features_10.rds")
)


cosine_dev_matrix_10 <-
  aml_data %>% 
  tof_feature_extract(
    method = "ddpr", 
    phenotypic_markers = c(all_of(c(SURFACE_MARKERS, TRX_FACTORS, "MPO"))), 
    signaling_markers = all_of(SIGNALING_MARKERS), 
    patient_var = patient, 
    condition_var = condition,
    cluster_var = cosine_cluster, 
    stimulation_var = stimulation, 
    stimulation_baseline = "Basal"
  )

write_rds(
  x = cosine_dev_matrix_10, 
  path = here::here("data", "cosine_dev_features_10.rds")
)
```

## CITRUS Extraction

``` r
mahalanobis_citrus_matrix <- 
  aml_data %>% 
  tof_feature_extract(
    method = "citrus", 
    phenotypic_markers = c(all_of(c(SURFACE_MARKERS, TRX_FACTORS, "MPO"))), 
    signaling_markers = all_of(SIGNALING_MARKERS), 
    patient_var = patient, 
    condition_var = condition,
    cluster_var = mahalanobis_cluster, 
    stimulation_var = stimulation, 
    stimulation_baseline = "Basal"
  )

write_rds(
  x = mahalanobis_citrus_matrix, 
  path = here::here("data", "mahalanobis_citrus_features.rds")
)

cosine_citrus_matrix <-
  aml_data %>% 
  tof_feature_extract(
    method = "citrus", 
    phenotypic_markers = c(all_of(c(SURFACE_MARKERS, TRX_FACTORS, "MPO"))), 
    signaling_markers = all_of(SIGNALING_MARKERS), 
    patient_var = patient, 
    condition_var = condition,
    cluster_var = cosine_cluster, 
    stimulation_var = stimulation, 
    stimulation_baseline = "Basal"
  )

write_rds(
  x = cosine_citrus_matrix, 
  path = here::here("data", "cosine_citrus_features.rds")
)
```

## Other Extraction
