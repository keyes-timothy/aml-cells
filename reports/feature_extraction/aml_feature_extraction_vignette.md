AML Feature Extraction Vignette
================
tkeyes
2021-02-19

  - [Executive Summary](#executive-summary)
  - [Read in data](#read-in-data)
  - [DDPR](#ddpr)
      - [Helper functions](#helper-functions)
      - [Calculate DDPR feature
        matrices](#calculate-ddpr-feature-matrices)
      - [Calculate abundance + central tendency feature
        matrices](#calculate-abundance-central-tendency-feature-matrices)
      - [Calculate “CITRUS” feature
        matrices](#calculate-citrus-feature-matrices)
      - [Calculate feature matrices (with EMD for
        signaling)](#calculate-feature-matrices-with-emd-for-signaling)
          - [Calculate feature matrices with JS
            Index](#calculate-feature-matrices-with-js-index)

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
    "emdist",
    "philentropy", 
    "tidytof"
  )

source(here::here("scripts", "setup", "aml_utils.R"))
call_libraries(libraries)
source(here::here("scripts", "feature_extraction", "feature_extraction_utils.R"))


# Parameters
set_global_variables(locale = "galaxia")
md_path <- here::here("data-raw", "AML_metadata.xlsx")

CLUSTER_OUTPUT <- here::here("data", "cancer_data_clustered.rds")

# set up aml marker and patient information
marker_setup()
patient_setup()

# Misc globals for running different parts of this vignette
is_sampled <- FALSE
aml_path <- here::here("data", "aml_data_clustered.rds")
all_path <- here::here("data", "all_data_clustered.rds")
sample_string <- ""

if (is_sampled) { 
  aml_path <- aml_path %>% 
    str_replace(".rds", "_sampled.rds")
  all_path <- all_path %>% 
    str_replace(".rds", "_sampled.rds")
  sample_string <- "_sampled"
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

If we want to use sampled data (for development), we take 10% of the
cells that were read in from each file.

# DDPR

The DDPR feature matrix has 3 types of features:

  - Proportion of cells in each cell population *expanded* relative to
    healthy controls.
  - Mean expression of lineage markers in each of these sub-populations
  - Proportion of cells that are “positive” for each of the measured
    signaling markers (which we implement as being above a certain
    threshold).

To get to each of these feature types, we can create some functions that
get us most of the way there and then wrap them in an overall
`tof_extract_features` function that allows us to specify ddpr as the
method we want to use.

## Helper functions

First, we source some helper functions that haven’t yet been
incorporated into `tidytof.`

``` r
source(here::here("scripts", "feature_extraction", "feature_extraction_utils.R"))
```

## Calculate DDPR feature matrices

The first step to compute the DDPR feature matrices is to separate the
cells that come from healthy samples and those that come from cancer
samples.

``` r
healthy_tibble <- 
  aml_data %>% 
  filter(condition == "healthy")


cancer_tibble <- 
  aml_data %>% 
  filter(condition != "healthy") %>% 
  group_by(patient, condition)
```

Then, we can find the DDPR feature matrices for each patient and
condition using a variety of signaling thresholds using the helper
function `aml_iterate_ddpr_extraction` in `aml_utils.R`.

We do this below using a variety of clustering techniques to define the
clusters from which we are extracting ddpr features including the
following:

  - developmental classification (using mahalanobis distance)
  - flowsom clustering (using surface markers)
  - flowsom clustering using the mahalanobis distances from single-cell
    classification
  - flowsom clustering using the cosine distances from single-cell
    classification

The result, a tibble with a list-col containing the extracted features,
will be saved to our data folder.

``` r
ddpr_feature_matrices <- 
  tibble(
    threshold = c(5, 10, 20, 50, 100), 
    features_mahalanobis = 
      map(
        .x = threshold, 
        .f = ~ 
          aml_iterate_ddpr_extraction(
            cluster_col = mahalanobis_cluster, 
            threshold = .x
          )
      ),
    features_flowsom_surface = 
      map(
        .x = threshold, 
        .f = ~ 
          aml_iterate_ddpr_extraction(
            cluster_col = cluster_flowSOM_surface, 
            threshold = .x
          )
      ),
    features_flowsom_mah = 
      map(
        .x = threshold, 
        .f = ~ 
          aml_iterate_ddpr_extraction(
            cluster_col = cluster_flowSOM_mah, 
            threshold = .x
          )
      ), 
    features_flowsom_cos = 
      map(
        .x = threshold, 
        .f = ~ 
          aml_iterate_ddpr_extraction(
            cluster_col = cluster_flowSOM_cos, 
            threshold = .x
          )
      )
  ) %>% 
  pivot_longer(
    cols = starts_with("features_"), 
    names_to = "cluster_type", 
    values_to = "features", 
    names_prefix = "features_"
  )



ddpr_feature_matrices
```

    ## # A tibble: 20 x 3
    ##    threshold cluster_type    features             
    ##        <dbl> <chr>           <list>               
    ##  1         5 mahalanobis     <tibble [42 × 605]>  
    ##  2         5 flowsom_surface <tibble [42 × 2,012]>
    ##  3         5 flowsom_mah     <tibble [42 × 2,012]>
    ##  4         5 flowsom_cos     <tibble [42 × 2,012]>
    ##  5        10 mahalanobis     <tibble [42 × 605]>  
    ##  6        10 flowsom_surface <tibble [42 × 2,012]>
    ##  7        10 flowsom_mah     <tibble [42 × 2,012]>
    ##  8        10 flowsom_cos     <tibble [42 × 2,012]>
    ##  9        20 mahalanobis     <tibble [42 × 605]>  
    ## 10        20 flowsom_surface <tibble [42 × 2,012]>
    ## 11        20 flowsom_mah     <tibble [42 × 2,012]>
    ## 12        20 flowsom_cos     <tibble [42 × 2,012]>
    ## 13        50 mahalanobis     <tibble [42 × 605]>  
    ## 14        50 flowsom_surface <tibble [42 × 2,012]>
    ## 15        50 flowsom_mah     <tibble [42 × 2,012]>
    ## 16        50 flowsom_cos     <tibble [42 × 2,012]>
    ## 17       100 mahalanobis     <tibble [42 × 605]>  
    ## 18       100 flowsom_surface <tibble [42 × 2,012]>
    ## 19       100 flowsom_mah     <tibble [42 × 2,012]>
    ## 20       100 flowsom_cos     <tibble [42 × 2,012]>

``` r
file_name <- str_c("ddpr_feature_matrices", sample_string, ".rds")

write_rds(
  x = ddpr_feature_matrices, 
  file = 
    here::here(
      "data", 
      "patient-level_feature_matrices", 
      file_name
    )
)
```

## Calculate abundance + central tendency feature matrices

We can also extract feature matrices from our single-cell data that
simply extract the mean of each of our markers in each clustered cell
type within each stimulation condition. For good measure, we can also
include the proportion of each cluster in these matrices.

For this, we will use the same 4 clustering strategies as above. First
we start with finding the means…

``` r
mean_feature_matrices <-
  tibble(
    features_mahalanobis = 
      list(
        aml_extract_lineage_features(
          cluster_col = mahalanobis_cluster, 
          central_tendency = mean
        )
      ),
    features_flowsom_surface = 
      list(
        aml_extract_lineage_features(
          cluster_col = cluster_flowSOM_surface, 
          central_tendency = mean
        )
      ),
    features_flowsom_mah = 
      list(
        aml_extract_lineage_features(
          cluster_col = cluster_flowSOM_mah, 
          central_tendency = mean
        )
      ), 
    features_flowsom_cos = 
      list(
        aml_extract_lineage_features(
          cluster_col = cluster_flowSOM_cos, 
          central_tendency = mean
        )
      )
  ) %>% 
  pivot_longer(
    cols = starts_with("features_"), 
    names_to = "cluster_type", 
    values_to = "features", 
    names_prefix = "features_"
  )
```

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

``` r
mean_feature_matrices
```

    ## # A tibble: 4 x 2
    ##   cluster_type    features             
    ##   <chr>           <list>               
    ## 1 mahalanobis     <tibble [46 × 371]>  
    ## 2 flowsom_surface <tibble [46 × 1,232]>
    ## 3 flowsom_mah     <tibble [46 × 1,232]>
    ## 4 flowsom_cos     <tibble [46 × 1,232]>

…and then we find the abundances.

``` r
abundance_feature_matrices <-
  tibble(
    features_mahalanobis = 
      list(
        extract_feature_proportion(
          tof_tibble = group_by(aml_data, patient, condition),
          cluster_col = mahalanobis_cluster
        )
      ),
    features_flowsom_surface = 
      list(
        extract_feature_proportion(
          tof_tibble = group_by(aml_data, patient, condition),
          cluster_col = cluster_flowSOM_surface 
        )
      ),
    features_flowsom_mah = 
      list(
        extract_feature_proportion(
          tof_tibble = group_by(aml_data, patient, condition),
          cluster_col = cluster_flowSOM_mah 
        )
      ), 
    features_flowsom_cos = 
      list(
        extract_feature_proportion(
          tof_tibble = group_by(aml_data, patient, condition),
          cluster_col = cluster_flowSOM_cos 
        )
      )
  ) %>% 
  pivot_longer(
    cols = starts_with("features_"), 
    names_to = "cluster_type", 
    values_to = "features", 
    names_prefix = "features_"
  ) %>%
  rename(abundances = features)
```

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

And then we combine the mean and abundance features:

``` r
final_features <- 
  left_join(
    mean_feature_matrices, 
    abundance_feature_matrices, 
    by = "cluster_type"
  ) %>% 
  transmute(
    cluster_type, 
    final_features = 
      map2(features, abundances, left_join, by = c("patient", "condition"))
  )


file_name <- str_c("mean_feature_matrices", sample_string, ".rds")

write_rds(
  x = final_features, 
  file = 
    here::here(
      "data", 
      "patient-level_feature_matrices", 
      file_name
    )
)
```

## Calculate “CITRUS” feature matrices

The citrus feature matrix is similar to the one taken directly above,
but takes into account which stimulation condition each cell is in (and
separates each stimulation condition into its own unique feature for
each cell population and marker).

For now, we won’t worry about calculating the CITRUS feature matrix
because with such a small number of samples, I am worried about
overfitting (with a very high number of features)

``` r
citrus_feature_matrices <- NULL
```

## Calculate feature matrices (with EMD for signaling)

Instead of calculating a measurement of central tendency for each
signaling feature (for each stimulation), can also look at measures of
how different the entire distributions are between the basal condition
and each of the stimulation conditions.

One way of doing this is to calculate the Earth Mover’s Distance(EMD)
between the distribution of cells in the basal condition and the
distribution of cells in each stimulation condition within a given
patient. This will yield a “score” measuring how different the
distributions (one score for each non-basal stimulation condition).

Note that calculating the EMD is very computationally expensive. A
single EMD calculation for one set of clusters (on heavily sampled data)
takes about 14 minutes locally.

``` r
emd_feature_matrices <- 
  list(
    ddpr_mahalanobis = 
      aml_extract_lineage_emd(mahalanobis_cluster, mean) %>% 
      ungroup(), 
    ddpr_cos = 
      aml_extract_lineage_emd(cosine_cluster, mean) %>% 
      ungroup(), 
    flowsom_surface = 
      aml_extract_lineage_emd(cluster_flowSOM_surface, mean) %>% 
      ungroup(), 
    flowsom_mah = 
      aml_extract_lineage_emd(cluster_flowSOM_mah, mean) %>% 
      ungroup(), 
    flowsom_cos = 
      aml_extract_lineage_emd(cluster_flowSOM_cos, mean) %>% 
      ungroup()
  )
```

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Joining, by = c("patient", "condition")

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Joining, by = c("patient", "condition")

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Joining, by = c("patient", "condition")

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Joining, by = c("patient", "condition")

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Joining, by = c("patient", "condition")

``` r
final_features <- 
  tibble(
    cluster_type = names(emd_feature_matrices), 
    features = emd_feature_matrices
  )

file_name <- str_c("emd_feature_matrices", sample_string, ".rds")


write_rds(
  x = final_features, 
  file = 
    here::here(
      "data", 
      "patient-level_feature_matrices", 
      file_name
    )
)
```

### Calculate feature matrices with JS Index

Another measure of the distance between two distributions is the
Jensen-Shannon index. We can calculate it for all of our samples similar
to how we calculated EMD…

The JSI is not quite as computationally expensive as the EMD.

``` r
js_feature_matrices <- 
  list(
    ddpr_mahalanobis = 
      aml_extract_lineage_js(mahalanobis_cluster, mean) %>%
      ungroup(),
    ddpr_cos =
      aml_extract_lineage_js(cosine_cluster, mean) %>%
      ungroup(),
    flowsom_surface = 
      aml_extract_lineage_js(cluster_flowSOM_surface, mean) %>% 
      ungroup(),
    flowsom_mah =
      aml_extract_lineage_js(cluster_flowSOM_mah, mean) %>%
      ungroup(),
    flowsom_cos =
      aml_extract_lineage_js(cluster_flowSOM_cos, mean) %>%
      ungroup()
  )
```

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## Joining, by = c("patient", "condition")

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## Joining, by = c("patient", "condition")

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## Joining, by = c("patient", "condition")

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## Joining, by = c("patient", "condition")

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## `summarise()` has grouped output by 'patient', 'condition'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `patient`, `condition`

    ## Joining, by = c("patient", "condition")

``` r
js_tibble <- 
  tibble(
    cluster_type = names(js_feature_matrices), 
    features = js_feature_matrices
  )

js_tibble
```

    ## # A tibble: 5 x 2
    ##   cluster_type     features             
    ##   <chr>            <named list>         
    ## 1 ddpr_mahalanobis <tibble [46 × 641]>  
    ## 2 ddpr_cos         <tibble [46 × 641]>  
    ## 3 flowsom_surface  <tibble [46 × 2,132]>
    ## 4 flowsom_mah      <tibble [46 × 2,132]>
    ## 5 flowsom_cos      <tibble [46 × 2,132]>

``` r
file_name <- str_c("js_feature_matrices", sample_string, ".rds")


write_rds(
  x = js_tibble,
  file =
    here::here(
      "data",
      "patient-level_feature_matrices",
      file_name
    )
)
```
