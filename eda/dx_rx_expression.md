Diagnosis vs. Relapse Expression Analysis
================
Timothy Keyes
2020-04-03

  - [Basic setup](#basic-setup)
  - [Cluster to cluster across
    patients](#cluster-to-cluster-across-patients)
      - [CD45](#cd45)
      - [CD34](#cd34)
      - [CD38](#cd38)
      - [CD90](#cd90)
      - [CD93](#cd93)
      - [CD99](#cd99)

``` r
# Libraries
library(tidyverse)

# Parameters
input_path <- here::here("data", "AML_matrix_clustered.rds")

CLASSIFIER_POPULATIONS <- 
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

#===============================================================================

#source necessary files
source('~/GitHub/aml-cells/scripts/setup/marker_setup.R')
source('~/GitHub/aml-cells/scripts/setup/patient_setup.R')
```

## Basic setup

Algorithm (in words):

  - Read in previously pre-processed data
      - arcsinh transform
      - compensation
      - filtered on `time` variable to remove weird CD34/CD38
        “debarcode-gate”
  - filter to only include paired samples (with both `Dx` and `Rx`
    present as well as healthy controls)

<!-- end list -->

``` r
#### set up global varibles
marker_setup()
patient_setup()

recode_vars <- 
  tibble(
    from = 1:length(CLASSIFIER_POPULATIONS), 
    to = CLASSIFIER_POPULATIONS
  ) %>% 
  deframe()

#read in paired data (and healthy data) 
aml_data <- 
  input_path %>% 
  read_rds() %>%
  mutate(
    patient = 
      patient %>% 
      as.character() %>% 
      str_to_lower() %>% 
      if_else(. == "pastrp", "pasrtp", .), 
    Mah.cluster = 
      recode(Mah.cluster, !!! recode_vars) %>% 
      factor(levels = CLASSIFIER_POPULATIONS)
  ) %>% 
  mutate_at(
    .vars = vars(phenograph.metacluster, FlowSOM.cluster, FlowSOM.metacluster),
    ~ (.) %>% 
      as.character() %>% 
      fct_reorder(., CD34, .fun = mean, .desc = TRUE)
  ) %>% 
  filter(patient %in% c(PAIRED_PATIENTS, HEALTHY_CONTROLS)) %>% 
  rename_at(
    vars(everything()), 
    ~ (.) %>% 
      str_to_lower() %>% 
      str_replace_all(pattern = "[:punct:]", replacement = "_")
  )
```

## Cluster to cluster across patients

``` r
cluster_column <- "mah_cluster"
marker <- "cd45"

clusters_across_patients <- function(cluster_column, marker) {
  
  map(
    levels(aml_data[[cluster_column]]), 
    function(cluster_type) {
      patient_levels <- 
        aml_data %>% 
        filter(!! sym(cluster_column) == cluster_type) %>% 
        group_by(
          metacondition = if_else(condition == "Healthy", "Healthy", "Cancer"),
          patient
        ) %>% 
        summarize(median = median(!! sym(marker))) %>% 
        ungroup() %>% 
        arrange(metacondition, median) %>% 
        pull(patient)
      
      aml_data %>% 
        filter(!! sym(cluster_column) %in% cluster_type) %>% 
        mutate(
          metacondition = 
            if_else(condition == "Healthy", "Healthy", "Cancer") %>% 
            as.factor(), 
          patient = factor(patient, levels = patient_levels)
        ) %>% 
        ggplot(aes(x = patient, y = !! sym(marker), fill = condition)) + 
        geom_violin(draw_quantiles = 0.5) + 
        labs(
          title = 
            str_c(
              str_to_upper(marker), 
              " expression in ", 
              cluster_type, 
              " cells across patients"
            ), 
          fill = NULL
        )
    }
  )
}
```

### CD45

``` r
cd45_mah_cluster_plots <- clusters_across_patients("mah_cluster", "cd45")

cd45_mah_cluster_plots %>% 
  walk(print)
```

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->![](dx_rx_expression_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->![](dx_rx_expression_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->![](dx_rx_expression_files/figure-gfm/unnamed-chunk-4-4.png)<!-- -->![](dx_rx_expression_files/figure-gfm/unnamed-chunk-4-5.png)<!-- -->![](dx_rx_expression_files/figure-gfm/unnamed-chunk-4-6.png)<!-- -->![](dx_rx_expression_files/figure-gfm/unnamed-chunk-4-7.png)<!-- -->![](dx_rx_expression_files/figure-gfm/unnamed-chunk-4-8.png)<!-- -->![](dx_rx_expression_files/figure-gfm/unnamed-chunk-4-9.png)<!-- -->

-----

### CD34

``` r
cd34_mah_cluster_plots <- clusters_across_patients("mah_cluster", "cd34")

cd34_mah_cluster_plots %>% 
  walk(print)
```

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->![](dx_rx_expression_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->![](dx_rx_expression_files/figure-gfm/unnamed-chunk-5-5.png)<!-- -->![](dx_rx_expression_files/figure-gfm/unnamed-chunk-5-6.png)<!-- -->![](dx_rx_expression_files/figure-gfm/unnamed-chunk-5-7.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-5-8.png)<!-- -->![](dx_rx_expression_files/figure-gfm/unnamed-chunk-5-9.png)<!-- -->

-----

### CD38

``` r
cd38_mah_cluster_plots <- clusters_across_patients("mah_cluster", "cd38")

cd38_mah_cluster_plots %>% 
  walk(print)
```

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->![](dx_rx_expression_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-6-4.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-6-5.png)<!-- -->![](dx_rx_expression_files/figure-gfm/unnamed-chunk-6-6.png)<!-- -->![](dx_rx_expression_files/figure-gfm/unnamed-chunk-6-7.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-6-8.png)<!-- -->![](dx_rx_expression_files/figure-gfm/unnamed-chunk-6-9.png)<!-- -->

-----

### CD90

``` r
cd90_mah_cluster_plots <- clusters_across_patients("mah_cluster", "cd90")

cd90_mah_cluster_plots %>% 
  walk(print)
```

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-7-3.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-7-4.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-7-5.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-7-6.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-7-7.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-7-8.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-7-9.png)<!-- -->

-----

### CD93

``` r
cd93_mah_cluster_plots <- clusters_across_patients("mah_cluster", "cd93")

cd93_mah_cluster_plots %>% 
  walk(print)
```

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-8-4.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-8-5.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-8-6.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-8-7.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-8-8.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-8-9.png)<!-- -->

-----

### CD99

``` r
cd99_mah_cluster_plots <- clusters_across_patients("mah_cluster", "cd99")

cd99_mah_cluster_plots %>% 
  walk(print)
```

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->![](dx_rx_expression_files/figure-gfm/unnamed-chunk-9-3.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values
    
    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-9-4.png)<!-- -->![](dx_rx_expression_files/figure-gfm/unnamed-chunk-9-5.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-9-6.png)<!-- -->![](dx_rx_expression_files/figure-gfm/unnamed-chunk-9-7.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-9-8.png)<!-- -->

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](dx_rx_expression_files/figure-gfm/unnamed-chunk-9-9.png)<!-- -->
