---
title: "Clustering AML vignette"
author: tkeyes
date: "2021-08-24"
output: 
  github_document:
    toc: true
---



# Executive Summary

In this vignette, we perform clustering on the AML and ALL datasets using several clustering 
algorithms and interrogate some of their basic properties. 


```r
# Libraries
libraries <- 
  c(
    "tidyverse", 
    "rlang", 
    "ggthemes", 
    "doParallel", 
    "ggrepel", 
    "FlowSOM", 
    "tidytof"
  )

# Parameters - may need to change on scg machines
source(here::here("scripts", "setup", "aml_utils.R"))
call_libraries(libraries)

set_global_variables(locale = "galaxia")

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

CLUSTER_OUTPUT <- 
  here::here("data", "aml_data_clustered.rds")

# set up aml marker and patient information
marker_setup()
patient_setup()

# globals for running different parts of this vignette
is_sampled <- FALSE

# paths 
healthy_myeloid_path <- here::here("data-raw", "healthy_myeloid")
healthy_B_path <- here::here("data-raw", "healthy_DDPR")

input_aml_path <- 
  here::here(
    "data", 
    if_else(is_sampled, "aml_data_downsampled_p.rds", "aml_data_p.rds")
  )
```

# Set-up


```r
# read in AML data - already processed and arcsinh transformed
aml_data <- 
  input_aml_path %>% 
  read_rds()


# load healthy myeloid and B-cell data
healthy_myeloid_data <- 
  healthy_myeloid_path %>% 
  tof_read_data(sep = "_") %>% 
  tof_preprocess()
```

```
## Error in tof_read_data(., sep = "_"): object 'tof_tibble' not found
```



```r
aml_data
```

```
## # A tibble: 22,527,299 × 52
##    file_name       CD45  CD61  CD99 empty_I127 CD45RA  CD93 CD3_CD19  CCR2 CD117
##    <chr>          <dbl> <dbl> <dbl>      <dbl>  <dbl> <dbl>    <dbl> <dbl> <dbl>
##  1 Plate1_BM5721…  4.84 0.199 0.390      0.199  3.75  0.569    0.390 0.199 0.199
##  2 Plate1_BM5721…  5.82 0.199 3.81       0.199  3.18  0.390    0.199 0.199 0.199
##  3 Plate1_BM5721…  5.16 3.64  3.69       0.199  1.94  2.89     0.569 4.20  0.199
##  4 Plate1_BM5721…  4.07 2.46  2.31       0.199  1.53  0.390    0.199 3.43  0.199
##  5 Plate1_BM5721…  4.61 0.199 3.42       0.199  4.45  0.199    0.199 0.199 0.199
##  6 Plate1_BM5721…  4.70 0.199 3.36       0.199  3.79  0.199    1.14  0.569 0.199
##  7 Plate1_BM5721…  5.75 1.44  3.00       1.44   2.19  2.56     0.390 3.85  0.199
##  8 Plate1_BM5721…  5.95 0.199 0.881      0.199  3.48  0.390    0.199 0.199 0.199
##  9 Plate1_BM5721…  5.12 0.199 1.14       0.199  2.94  0.199    1.61  0.199 0.199
## 10 Plate1_BM5721…  5.13 0.199 2.87       0.199  0.733 0.733    0.569 4.11  0.199
## # … with 22,527,289 more rows, and 42 more variables: CD123 <dbl>, CD64 <dbl>,
## #   CD90 <dbl>, CD38 <dbl>, CD34 <dbl>, CEBPa <dbl>, pSTAT5 <dbl>, CD11c <dbl>,
## #   CD13 <dbl>, pAkt <dbl>, TIM-3 <dbl>, CD56 <dbl>, CD10 <dbl>, PU.1 <dbl>,
## #   CD33 <dbl>, PD-L1 <dbl>, CD14 <dbl>, caspase-3 <dbl>, GATA-1 <dbl>,
## #   pSTAT3 <dbl>, CD41 <dbl>, CD16 <dbl>, CD68 <dbl>, MPO <dbl>, pErk <dbl>,
## #   CD47 <dbl>, CD135 <dbl>, CD109 <dbl>, pS6 <dbl>, CD49f <dbl>, HLA-DR <dbl>,
## #   CD71 <dbl>, pCreb <dbl>, empty_Ir191 <dbl>, empty_Ir193 <dbl>, …
```

```r
healthy_myeloid_data
```

```
## Error in eval(expr, envir, enclos): object 'healthy_myeloid_data' not found
```


We have to some quick tidying of the healthy data channel names: 


```r
my_channels <- 
  healthy_myeloid_data %>% 
  select(contains("_"), -Event_length_length, -`cisplatin_195Pt`) %>% 
  colnames() %>% 
  if_else(
    str_detect(., "empty") | !str_detect(., "_"), 
    ., 
    str_extract(., pattern = ".+_") %>% 
      str_sub(end = -2L)
  )
```

```
## Error in select(., contains("_"), -Event_length_length, -cisplatin_195Pt): object 'healthy_myeloid_data' not found
```

```r
my_channels <- 
  my_channels[my_channels %in% colnames(aml_data)]
```

```
## Error in eval(expr, envir, enclos): object 'my_channels' not found
```

```r
healthy_myeloid_data <- 
  healthy_myeloid_data %>% 
  rename_with(
    .fn =
      function(x) {
        if_else(
          str_detect(x, "empty") | !str_detect(x, "_"),
          x, 
          str_extract(x, pattern = ".+_") %>% 
            str_sub(end = -2L)
        )
      }
  ) %>% 
  select(file, all_of(my_channels)) %>% 
  rename(file_name = file)
```

```
## Error in rename_with(., .fn = function(x) {: object 'healthy_myeloid_data' not found
```


Then we have to extract some metadata from the file name for each classifier population


```r
healthy_myeloid_data <- 
  healthy_myeloid_data %>% 
  mutate(
    population = 
      file_name %>% 
      str_remove(pattern = ".+_") %>% 
      str_remove(pattern = ".fcs")
  )
```

```
## Error in mutate(., population = file_name %>% str_remove(pattern = ".+_") %>% : object 'healthy_myeloid_data' not found
```

```r
healthy_myeloid_data
```

```
## Error in eval(expr, envir, enclos): object 'healthy_myeloid_data' not found
```


# Developmental Classifier

## Build the classifier


```r
# mahalanobis clusters 
mahalanobis_clusters <-
  healthy_myeloid_data %>% 
  tof_cluster_ddpr(
    cancer_tibble = aml_data, 
    healthy_cell_labels = healthy_myeloid_data$population, 
    cluster_cols = all_of(CLASSIFIER_MARKERS),
    distance_function = "mahalanobis", 
    num_cores = 12L, 
    parallel_cols = patient,
    return_distances = TRUE
  )
```

```
## Error in dplyr::mutate(., population = healthy_cell_labels): object 'healthy_myeloid_data' not found
```

```r
mahalanobis_clusters
```

```
## Error in eval(expr, envir, enclos): object 'mahalanobis_clusters' not found
```

```r
# calculate cluster-wise accuracy in healthies
mahalanobis_clusters_healthy <-
  healthy_myeloid_data %>% 
  tof_cluster_ddpr(
    cancer_tibble = healthy_myeloid_data, 
    healthy_cell_labels = healthy_myeloid_data$population, 
    cluster_cols = all_of(CLASSIFIER_MARKERS),
    distance_function = "mahalanobis", 
    num_cores = 6L, 
    parallel_cols = file_name,
    return_distances = TRUE
  )
```

```
## Error in dplyr::mutate(., population = healthy_cell_labels): object 'healthy_myeloid_data' not found
```

```r
healthy_myeloid_data %>% 
  bind_cols(mahalanobis_clusters_healthy) %>% 
  group_by(population) %>% 
  summarize(prop_correct = sum(population == mahalanobis_cluster) / n()) %>% 
  arrange(-prop_correct)
```

```
## Error in list2(...): object 'healthy_myeloid_data' not found
```



```r
# cosine clusters
cosine_clusters <-
  healthy_myeloid_data %>% 
  tof_cluster_ddpr(
    cancer_tibble = aml_data, 
    healthy_cell_labels = healthy_myeloid_data$population, 
    cluster_cols = all_of(CLASSIFIER_MARKERS),
    distance_function = "cosine", 
    num_cores = 12L, 
    parallel_cols = patient,
    return_distances = TRUE
  )
```

```
## Error in dplyr::mutate(., population = healthy_cell_labels): object 'healthy_myeloid_data' not found
```

```r
cosine_clusters
```

```
## Error in eval(expr, envir, enclos): object 'cosine_clusters' not found
```

```r
# calculate cluster-wise accuracy in healthies
cosine_clusters_healthy <-
  healthy_myeloid_data %>% 
  tof_cluster_ddpr(
    cancer_tibble = healthy_myeloid_data, 
    healthy_cell_labels = healthy_myeloid_data$population, 
    cluster_cols = all_of(CLASSIFIER_MARKERS),
    distance_function = "cosine", 
    num_cores = 6L, 
    parallel_cols = file_name,
    return_distances = TRUE
  )
```

```
## Error in dplyr::mutate(., population = healthy_cell_labels): object 'healthy_myeloid_data' not found
```

```r
healthy_myeloid_data %>% 
  bind_cols(cosine_clusters_healthy) %>% 
  group_by(population) %>% 
  summarize(prop_correct = sum(population == cosine_cluster) / n()) %>% 
  arrange(-prop_correct)
```

```
## Error in list2(...): object 'healthy_myeloid_data' not found
```

pearson clusters


```r
# cosine clusters
pearson_clusters <-
  healthy_myeloid_data %>% 
  tof_cluster_ddpr(
    cancer_tibble = aml_data, 
    healthy_cell_labels = healthy_myeloid_data$population, 
    cluster_cols = all_of(CLASSIFIER_MARKERS),
    distance_function = "pearson", 
    num_cores = 12L, 
    parallel_cols = patient,
    return_distances = TRUE
  )
```

```
## Error in dplyr::mutate(., population = healthy_cell_labels): object 'healthy_myeloid_data' not found
```

```r
pearson_clusters
```

```
## Error in eval(expr, envir, enclos): object 'pearson_clusters' not found
```

```r
# calculate cluster-wise accuracy in healthies
pearson_clusters_healthy <-
  healthy_myeloid_data %>% 
  tof_cluster_ddpr(
    cancer_tibble = healthy_myeloid_data, 
    healthy_cell_labels = healthy_myeloid_data$population, 
    cluster_cols = all_of(CLASSIFIER_MARKERS),
    distance_function = "pearson", 
    num_cores = 6L, 
    parallel_cols = file_name,
    return_distances = TRUE
  )
```

```
## Error in dplyr::mutate(., population = healthy_cell_labels): object 'healthy_myeloid_data' not found
```

```r
healthy_myeloid_data %>% 
  bind_cols(pearson_clusters_healthy) %>% 
  group_by(population) %>% 
  summarize(prop_correct = sum(population == pearson_cluster) / n()) %>% 
  arrange(-prop_correct)
```

```
## Error in list2(...): object 'healthy_myeloid_data' not found
```

# Phenograph


```r
aml_data_nested <- 
  aml_data %>% 
  group_by(patient) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(num_cells = map_int(.x = data, .f = nrow))

head(aml_data_nested)
```

```
## # A tibble: 6 × 3
##   patient data                      num_cells
##   <chr>   <list>                        <int>
## 1 bm5721  <tibble [298,101 × 51]>      298101
## 2 papwhs  <tibble [428,579 × 51]>      428579
## 3 parbiu  <tibble [2,397,161 × 51]>   2397161
## 4 bm6631  <tibble [233,606 × 51]>      233606
## 5 pasgwh  <tibble [800,415 × 51]>      800415
## 6 pasxvc  <tibble [1,965,794 × 51]>   1965794
```


```r
temp <- 
  aml_data_nested %>% 
  slice_head(n = 2L) %>% 
  mutate(
    phenograph_clusters = 
      map(
        .x = data, 
        .f = tof_cluster_phenograph, 
        cluster_cols = any_of(SURFACE_MARKERS)
      )
  )
```

```
## Loading required namespace: Rphenograph
```

```
## Error: Problem with `mutate()` column `phenograph_clusters`.
## ℹ `phenograph_clusters = map(.x = data, .f = tof_cluster_phenograph, cluster_cols = any_of(SURFACE_MARKERS))`.
## ✖ This function requires the {Rphenograph} package. Install it with this code:
## 
##            if(!require(devtools)){
##               install.packages("devtools")
##            }
##            devtools::install_github("JinmiaoChenLab/Rphenograph")
```

```r
temp %>% 
  mutate(num_cells = map_int(data, nrow)) %>%
  select(-num_cells) %>% 
  unnest(cols = c(data, phenograph_clusters)) %>% 
  mutate(phenograph_cluster = str_c(patient, "_", phenograph_cluster)) %>% 
  arrange(as.numeric(cell_id)) %>% 
  count(phenograph_cluster)
```

```
## Error in mutate(., num_cells = map_int(data, nrow)): object 'temp' not found
```



```r
start <- Sys.time()
aml_data %>% 
  slice_sample(n = 1000) %>% 
  tof_cluster_phenograph(cluster_cols = all_of(SURFACE_MARKERS))
print(Sys.time() - start)

start <- Sys.time()
aml_data %>% 
  slice_sample(n = 5000) %>% 
  tof_cluster_phenograph(cluster_cols = all_of(SURFACE_MARKERS))
print(Sys.time() - start)

start <- Sys.time()
aml_data %>% 
  slice_sample(n = 10000) %>% 
  tof_cluster_phenograph(cluster_cols = all_of(SURFACE_MARKERS)) 
print(Sys.time() - start)

start <- Sys.time()
aml_data %>% 
  slice_sample(n = 20000) %>% 
  tof_cluster_phenograph(cluster_cols = all_of(SURFACE_MARKERS))
print(Sys.time() - start)

start <- Sys.time()
aml_data %>% 
  slice_sample(n = 40000) %>% 
  tof_cluster_phenograph(cluster_cols = all_of(SURFACE_MARKERS))
print(Sys.time() - start)

start <- Sys.time()
aml_data %>% 
  slice_sample(n = 80000) %>% 
  tof_cluster_phenograph(cluster_cols = all_of(SURFACE_MARKERS))
print(Sys.time() - start)
```


# FlowSOM


```r
# clustering with surface markers 
flowsom_cluster_surface <- 
  aml_data %>% 
  tof_cluster_flowsom(
    cluster_cols = any_of(SURFACE_MARKERS), 
    perform_metaclustering = FALSE
  )

flowsom_metacluster_surface <- 
  aml_data %>% 
  tof_cluster_flowsom(cluster_cols = any_of(SURFACE_MARKERS))

# clustering using mahalanobis distances as input to flowsom
flowsom_mahalanobis_clusters <- 
  mahalanobis_clusters %>%
  select(-mahalanobis_cluster) %>% 
  tof_cluster_flowsom(cluster_cols = everything())
```

```
## Error in select(., -mahalanobis_cluster): object 'mahalanobis_clusters' not found
```

```r
# clustering with cosine distances
flowsom_cosine_clusters <- 
  cosine_clusters %>%
  select(-cosine_cluster) %>% 
  tof_cluster_flowsom(cluster_cols = everything())
```

```
## Error in select(., -cosine_cluster): object 'cosine_clusters' not found
```

```r
# clustering with signaling markers
flowsom_signaling_clusters <- 
  aml_data %>% 
  tof_cluster_flowsom(cluster_cols = any_of(SIGNALING_MARKERS))

# clustering with all markers
flowsom_all_marker_clusters <- 
  aml_data %>% 
  tof_cluster_flowsom(cluster_cols = any_of(ALL_MARKERS))
```

Adding clusters to master feature matrices and saving

```r
aml_data_clustered <- 
  aml_data %>% 
  mutate(
    flowsom_cluster_surface = flowsom_cluster_surface, 
    flowsom_metacluster_surface = flowsom_metacluster_surface, 
    flowsom_cluster_signaling = flowsom_signaling_clusters, 
    flowsom_cluster_mahalanobis = flowsom_mahalanobis_clusters, 
    flowsom_cluster_cosine = flowsom_cosine_clusters, 
    flowsom_cluster_all_markers = flowsom_all_marker_clusters
  ) %>% 
  bind_cols(mahalanobis_clusters) %>% 
  bind_cols(cosine_clusters) %>% 
  bind_cols(pearson_clusters)
```

```
## Error: Problem with `mutate()` column `flowsom_cluster_mahalanobis`.
## ℹ `flowsom_cluster_mahalanobis = flowsom_mahalanobis_clusters`.
## ✖ object 'flowsom_mahalanobis_clusters' not found
```

```r
write_rds(
  x = aml_data_clustered, 
  file = 
    CLUSTER_OUTPUT %>% 
    if_else(
      is_sampled, 
      true = str_replace(string = ., pattern = ".rds", replacement = "_sampled.rds"), 
      false = .
    )
)
```

```
## Error in saveRDS(x, con, version = version, refhook = refhook): object 'aml_data_clustered' not found
```


