Clustering AML vignette
================
tkeyes
2020-09-15

  - [Executive Summary](#executive-summary)
  - [Set-up](#set-up)
  - [Developmental Classifier](#developmental-classifier)
      - [Build the classifiers](#build-the-classifiers)
      - [Apply the classifiers](#apply-the-classifiers)
  - [Phenograph](#phenograph)
  - [FlowSOM](#flowsom)
  - [Other kinds of clusters?](#other-kinds-of-clusters)

# Executive Summary

In this vignette, we perform clustering on the AML and ALL datasets
using several clustering algorithms and interrogate some of their basic
properties.

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
    "DataExplorer", 
    "ComplexHeatmap", 
    "foreach", 
    "doParallel", 
    "ggrepel", 
    "ggiraph", 
    "ggiraphExtra", 
    "FlowSOM"
  )

# Parameters
source(here::here("scripts", "setup", "aml_utils.R")) #may need to change on scg machines
call_libraries(libraries)

set_global_variables(locale = "galaxia")

tidyTOF_directory <- file.path("~", "GitHub", "tidyTOF")
md_path <- here::here("data-raw", "AML_metadata.xlsx")

ALL_CLASSIFIER_MARKERS <- 
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

CLUSTER_OUTPUT <- here::here("data", "cancer_data_clustered.rds")

# Sourcing tidyTOF functions
source_tidyTOF(tidyTOF_directory)

# set up aml marker and patient information
marker_setup()
patient_setup()

# globals for running different parts of this vignette
is_sampled <- TRUE

healthy_myeloid_path <- here::here("data-raw", "healthy_myeloid")
healthy_B_path <- here::here("data-raw", "healthy_DDPR")
```

# Set-up

``` r
# read in AML data - already processed and arcsinh transformed
aml_data <- load_cancer_data(cancer = "aml", sample = is_sampled)

#read in ALL data - already preprocessed and arcsinh transformed
all_data <- 
  load_cancer_data(cancer = "all", sample = is_sampled) %>% 
  rename(kappa_lambda = `IgL kappa`)

# load healthy myeloid and B-cell data
healthy_myeloid_data <- 
  healthy_myeloid_path %>% 
  tof_read_fcs() %>% 
  tof_preprocess()

healthy_B_data <- 
  healthy_B_path %>% 
  tof_read_fcs() %>% 
  tof_preprocess()
```

We have to some quick tidying of the healthy data channel names:

``` r
my_channels <- 
  healthy_myeloid_data %>% 
  select(contains("_"), -Event_length, -`195Pt_cisplatin`) %>% 
  colnames() %>% 
  str_replace(pattern = ".+_", replacement = "")

healthy_myeloid_data <- 
  healthy_myeloid_data %>% 
  rename_with(.fn = function(x) str_replace(x, pattern = ".+_", replacement = "")) %>% 
  select(all_of(my_channels)) %>% 
  rename(file_names = names, CD3_CD19 = CD19)

healthy_B_data <- 
  healthy_B_data %>% 
  rename(Ikaros = tIkaros, file_name = file_names, Tdt = TdT, kappa_lambda = Kappa_lambda)
```

Then we have to extract some metadata from the file name for each
classifier population

``` r
healthy_myeloid_data <- 
  healthy_myeloid_data %>% 
  mutate(
    population = 
      file_names %>% 
      str_remove(pattern = ".+_") %>% 
      str_remove(pattern = ".fcs")
  )

healthy_B_data <- 
  healthy_B_data %>% 
  mutate(
    population = 
      file_name %>% 
      str_remove("concat.+[:digit:]+\\.") %>% 
      str_remove(".fcs$")
  )
```

# Developmental Classifier

## Build the classifiers

``` r
# aml classifier 
aml_dev_classifier <-
  healthy_myeloid_data %>% 
  tof_classifier_build(
  population_vector = healthy_myeloid_data$population, 
  classifier_markers = CLASSIFIER_MARKERS
)
```

    ## Building classifier from healthy cells...

    ## Reshaping healthy data...

    ## Calculating mean and covariance matrix for all healthy populations

    ## Done! Returning classifier_fit object

``` r
all_dev_classifier <- 
  healthy_B_data %>% 
  tof_classifier_build(
    population_vector = healthy_B_data$population, 
    classifier_markers = ALL_CLASSIFIER_MARKERS
  )
```

    ## Building classifier from healthy cells...

    ## Reshaping healthy data...

    ## Calculating mean and covariance matrix for all healthy populations

    ## Done! Returning classifier_fit object

## Apply the classifiers

``` r
aml_dev_mah <- 
  aml_data %>% 
  tof_classifier_apply(
    classifier_fit = aml_dev_classifier, 
    num_cores = 10, 
    parallel_var = patient,
    dist_fun = "mahalanobis"
  )

all_dev_mah <- 
  all_data %>% 
  tof_classifier_apply(
    classifier_fit = all_dev_classifier, 
    num_cores = 10, 
    parallel_var = file_name,
    dist_fun = "mahalanobis"
  )

aml_dev_cos <- 
  aml_data %>% 
  tof_classifier_apply(
    classifier_fit = aml_dev_classifier, 
    num_cores = 10, 
    parallel_var = patient,
    dist_fun = "cosine"
  )

all_dev_cos <- 
  all_data %>% 
  tof_classifier_apply(
    classifier_fit = all_dev_classifier, 
    num_cores = 10, 
    parallel_var = file_name,
    dist_fun = "cosine"
  )
```

# Phenograph

# FlowSOM

``` r
# clustering with surface markers 
aml_flowSOM_surface <- 
  aml_data %>% 
  tof_cluster_flowSOM(clustering_markers = SURFACE_MARKERS, num_clusters = 30)
```

    ## Building SOM

    ## Mapping data to SOM

    ## Building MST

``` r
# clustering with mahalanobis distances 
mah_cols <- 
  aml_dev_mah %>% 
  select(-mahalanobis_cluster) %>% 
  colnames()

aml_flowSOM_mah <- 
  aml_dev_mah %>% 
  tof_cluster_flowSOM(clustering_markers = mah_cols, num_clusters = 30)
```

    ## Building SOM

    ## Mapping data to SOM

    ## Building MST

``` r
# clustering with cosine distances
cos_cols <- 
  aml_dev_cos %>% 
  select(-cosine_cluster) %>% 
  colnames()

aml_flowSOM_cos <- 
  aml_dev_cos %>% 
  tof_cluster_flowSOM(clustering_markers = cos_cols, num_clusters = 30)
```

    ## Building SOM

    ## Mapping data to SOM

    ## Building MST

``` r
# clustering with signaling markers
aml_flowSOM_signaling <- 
  aml_data %>% 
  tof_cluster_flowSOM(clustering_markers = SIGNALING_MARKERS, num_clusters = 30)
```

    ## Building SOM

    ## Mapping data to SOM

    ## Building MST

``` r
# clustering with all markers
aml_flowSOM_all_markers <- 
  aml_data %>% 
  tof_cluster_flowSOM(clustering_markers = ALL_MARKERS, num_clusters = 30)
```

    ## Building SOM

    ## Mapping data to SOM

    ## Building MST

``` r
# clustering with all features
```

Adding clusters to master feature matrices and saving

``` r
#aml

aml_data_clustered <- 
  aml_data %>% 
  mutate(
    cluster_flowSOM_surface = aml_flowSOM_surface, 
    cluster_flowSOM_signaling = aml_flowSOM_signaling, 
    cluster_flowSOM_mah = aml_flowSOM_mah, 
    cluster_flowSOM_cos = aml_flowSOM_cos, 
    cluster_flowSOM_all_markers = aml_flowSOM_all_markers
  ) %>% 
  bind_cols(aml_dev_mah) %>% 
  bind_cols(aml_dev_cos)

saveRDS(
  aml_data_clustered, 
  file = 
    CLUSTER_OUTPUT %>% 
    str_replace("cancer", "aml") %>% 
    if_else(
      is_sampled, 
      true = str_replace(string = ., pattern = ".rds", replacement = "_sampled.rds"), 
      false = .
    )
)

#all

all_data_clustered <- 
  all_data %>% 
  bind_cols(all_dev_mah) %>% 
  bind_cols(all_dev_cos)

saveRDS(
  all_data_clustered, 
  file = 
    CLUSTER_OUTPUT %>% 
    str_replace("cancer", "all") %>% 
    if_else(
      is_sampled, 
      true = str_replace(string = ., pattern = ".rds", replacement = "_sampled.rds"), 
      false = .
    )
)
```

# Other kinds of clusters?

Future directions: other clustering approaches…

  - K-means
  - Phenograph (slow\!)
  - X-shift (will need to implement de novo)
