Disease-specific CyTOF Analysis
================
tkeyes
2021-07-12

-   [Background](#background)
-   [Read in data](#read-in-data)
-   [Fit linear models (everything
    together)](#fit-linear-models-everything-together)
-   [build all-cell model with only classifier
    markers](#build-all-cell-model-with-only-classifier-markers)
-   [build all-cell model with only signaling
    markers](#build-all-cell-model-with-only-signaling-markers)
-   [build a separate linear model for each cell
    type](#build-a-separate-linear-model-for-each-cell-type)

``` r
# Libraries
libraries <- 
  c(
    "tidyverse", 
    "tidymodels", 
    "rlang", 
    "ggthemes", 
    "doParallel", 
    "vip", 
    "tidytof"
  )

source(here::here("scripts", "setup", "aml_utils.R"))
source(here::here("scripts", "setup", "disease_specific_utils.R")) 

call_libraries(libraries)

# Parameters
set_global_variables(locale = "galaxia")
input_path <- here::here("data", "aml_data_clustered_sampled.rds")

# set up aml marker and patient information
marker_setup()
patient_setup()
```

# Background

# Read in data

``` r
aml_data <- 
  input_path %>% 
  read_rds()

aml_data
```

    ## # A tibble: 2,252,598 x 76
    ##     CD45  CD61  CD99 `127I` CD45RA  CD93 CD3_CD19  CCR2 CD117 CD123  CD64  CD90
    ##    <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>    <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    ##  1  4.86 1.61  2.73   0.733   2.19 0.390    1.61  3.64  0.199 0.881 3.13  0.199
    ##  2  5.82 2.42  3.33   0.199   2.62 2.64     0.569 3.96  0.199 0.199 2.89  0.199
    ##  3  5.70 1.88  1.82   0.199   2.52 2.67     0.733 3.13  0.199 0.390 2.89  0.733
    ##  4  4.64 0.390 2.09   0.199   2.19 0.390    0.733 3.15  0.569 3.27  0.569 0.569
    ##  5  6.05 1.44  2.92   0.199   1.75 1.02     0.199 4.01  0.199 0.390 3.62  0.199
    ##  6  5.68 2.39  2.59   0.199   2.31 2.49     0.199 2.82  0.199 0.199 2.75  0.390
    ##  7  5.09 0.199 0.733  0.199   1.14 0.199    0.881 4.09  0.199 4.05  0.199 0.733
    ##  8  4.91 0.199 3.63   0.390   1.53 1.25     1.35  4.48  0.199 0.199 1.61  0.199
    ##  9  5.62 0.199 2.64   0.199   1.61 0.199    0.199 0.390 0.199 0.199 0.199 0.199
    ## 10  3.49 0.569 3.15   0.199   1.82 0.199    0.199 0.199 0.199 0.569 0.199 0.199
    ## # … with 2,252,588 more rows, and 64 more variables: CD38 <dbl>, CD34 <dbl>,
    ## #   CEBPa <dbl>, pSTAT5 <dbl>, CD11c <dbl>, CD13 <dbl>, pAkt <dbl>,
    ## #   TIM-3 <dbl>, CD56 <dbl>, CD10 <dbl>, PU.1 <dbl>, CD33 <dbl>, CD14 <dbl>,
    ## #   caspase-3 <dbl>, GATA-1 <dbl>, pSTAT3 <dbl>, CD41 <dbl>, CD16 <dbl>,
    ## #   CD68 <dbl>, MPO <dbl>, pErk <dbl>, CD47 <dbl>, CD135 <dbl>, CD109 <dbl>,
    ## #   pS6 <dbl>, CD49f <dbl>, HLA-DR <dbl>, CD71 <dbl>, pCreb <dbl>, 191Ir <dbl>,
    ## #   193Ir <dbl>, cisplatin <dbl>, CD11b <dbl>, file_names <chr>, plate <fct>,
    ## #   patient <chr>, stimulation <chr>, condition <chr>, cell_id <chr>,
    ## #   cluster_flowSOM_surface <chr>, cluster_flowSOM_signaling <chr>,
    ## #   cluster_flowSOM_mah <chr>, cluster_flowSOM_cos <chr>,
    ## #   cluster_flowSOM_all_markers <chr>, mahalanobis_HSC <dbl>,
    ## #   mahalanobis_MPP <dbl>, mahalanobis_CMP <dbl>, mahalanobis_GMP <dbl>,
    ## #   mahalanobis_MEP <dbl>, mahalanobis_Monocyte <dbl>, mahalanobis_DC <dbl>,
    ## #   mahalanobis_Macrophage <dbl>, mahalanobis_Thrombocyte <dbl>,
    ## #   mahalanobis_cluster <chr>, cosine_HSC <dbl[,1]>, cosine_MPP <dbl[,1]>,
    ## #   cosine_CMP <dbl[,1]>, cosine_GMP <dbl[,1]>, cosine_MEP <dbl[,1]>,
    ## #   cosine_Monocyte <dbl[,1]>, cosine_DC <dbl[,1]>,
    ## #   cosine_Macrophage <dbl[,1]>, cosine_Thrombocyte <dbl[,1]>,
    ## #   cosine_cluster <chr>

``` r
colnames(aml_data)
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
aml_sampled <- 
  aml_data %>% 
  group_by(patient) %>% 
  slice_sample(prop = 0.01) %>% 
  ungroup()

aml_sampled
```

    ## # A tibble: 22,508 x 76
    ##     CD45  CD61  CD99 `127I` CD45RA  CD93 CD3_CD19  CCR2 CD117 CD123  CD64  CD90
    ##    <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>    <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    ##  1  5.67 1.61  3.11   0.199   1.53 1.82     0.390 3.50  0.199 0.199 2.98  0.199
    ##  2  5.53 0.199 0.199  0.199   2.59 0.733    0.199 0.199 0.199 1.75  0.199 0.199
    ##  3  5.33 3.00  3.58   0.199   2.52 1.94     1.44  4.13  0.199 1.14  3.33  0.199
    ##  4  6.11 2.56  2.09   0.199   1.02 0.390    0.390 1.25  0.199 1.02  1.99  0.199
    ##  5  6.02 0.199 2.27   0.199   2.14 0.199    0.569 0.390 0.199 0.199 0.199 0.390
    ##  6  5.86 1.75  3.49   0.199   2.39 0.733    0.881 3.84  0.199 0.199 4.20  0.733
    ##  7  5.24 3.36  2.98   0.199   2.85 1.88     1.75  3.73  0.199 0.199 2.80  0.199
    ##  8  5.84 3.60  1.88   0.199   1.88 1.35     1.35  2.92  0.199 0.390 3.18  0.199
    ##  9  4.49 0.733 0.881  0.199   1.68 0.199    0.199 3.73  0.199 4.03  0.199 0.199
    ## 10  5.77 2.42  1.53   0.199   1.94 1.14     0.569 1.35  0.390 1.14  2.23  0.199
    ## # … with 22,498 more rows, and 64 more variables: CD38 <dbl>, CD34 <dbl>,
    ## #   CEBPa <dbl>, pSTAT5 <dbl>, CD11c <dbl>, CD13 <dbl>, pAkt <dbl>,
    ## #   TIM-3 <dbl>, CD56 <dbl>, CD10 <dbl>, PU.1 <dbl>, CD33 <dbl>, CD14 <dbl>,
    ## #   caspase-3 <dbl>, GATA-1 <dbl>, pSTAT3 <dbl>, CD41 <dbl>, CD16 <dbl>,
    ## #   CD68 <dbl>, MPO <dbl>, pErk <dbl>, CD47 <dbl>, CD135 <dbl>, CD109 <dbl>,
    ## #   pS6 <dbl>, CD49f <dbl>, HLA-DR <dbl>, CD71 <dbl>, pCreb <dbl>, 191Ir <dbl>,
    ## #   193Ir <dbl>, cisplatin <dbl>, CD11b <dbl>, file_names <chr>, plate <fct>,
    ## #   patient <chr>, stimulation <chr>, condition <chr>, cell_id <chr>,
    ## #   cluster_flowSOM_surface <chr>, cluster_flowSOM_signaling <chr>,
    ## #   cluster_flowSOM_mah <chr>, cluster_flowSOM_cos <chr>,
    ## #   cluster_flowSOM_all_markers <chr>, mahalanobis_HSC <dbl>,
    ## #   mahalanobis_MPP <dbl>, mahalanobis_CMP <dbl>, mahalanobis_GMP <dbl>,
    ## #   mahalanobis_MEP <dbl>, mahalanobis_Monocyte <dbl>, mahalanobis_DC <dbl>,
    ## #   mahalanobis_Macrophage <dbl>, mahalanobis_Thrombocyte <dbl>,
    ## #   mahalanobis_cluster <chr>, cosine_HSC <dbl[,1]>, cosine_MPP <dbl[,1]>,
    ## #   cosine_CMP <dbl[,1]>, cosine_GMP <dbl[,1]>, cosine_MEP <dbl[,1]>,
    ## #   cosine_Monocyte <dbl[,1]>, cosine_DC <dbl[,1]>,
    ## #   cosine_Macrophage <dbl[,1]>, cosine_Thrombocyte <dbl[,1]>,
    ## #   cosine_cluster <chr>

# Fit linear models (everything together)

fit all cells to the healthy subspace

``` r
# find healthies to be used as the subspace
healthy_sampled <- 
  aml_data %>% 
  filter(condition == "healthy") %>% 
  select(any_of(ALL_MARKERS)) %>% 
  select(-CD3_CD19)

healthy_colnames <- colnames(healthy_sampled) 

healthy_sampled <- 
  healthy_sampled %>% 
  t() %>% 
  as_tibble()
```

    ## Warning: The `x` argument of `as_tibble.matrix()` must have unique column names if `.name_repair` is omitted as of tibble 2.0.0.
    ## Using compatibility `.name_repair`.

``` r
healthy_sampled
```

    ## # A tibble: 40 x 136,272
    ##       V1    V2    V3    V4    V5    V6    V7    V8    V9   V10   V11   V12   V13
    ##    <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    ##  1 0.199 0.390 0.199 0.199 0.199 0.569 0.199 0.199 0.199 0.733 0.199 0.199 0.199
    ##  2 1.14  3.94  4.11  2.31  1.88  3.97  2.19  0.733 1.44  0.390 3.86  2.96  2.42 
    ##  3 3.13  4.19  4.59  1.25  4.53  4.17  2.05  4.03  4.23  0.199 4.73  4.07  2.85 
    ##  4 2.78  3.78  5.19  4.00  4.66  3.76  3.65  4.09  0.199 0.881 4.70  4.62  3.04 
    ##  5 2.31  4.26  4.99  0.199 4.35  4.07  0.199 0.881 0.199 0.199 4.39  3.44  1.53 
    ##  6 0.199 0.199 0.199 0.199 0.199 0.199 0.199 0.199 0.199 0.199 0.199 0.199 0.199
    ##  7 3.23  4.39  3.78  0.569 4.30  3.73  0.881 3.06  0.199 0.199 4.57  3.62  3.00 
    ##  8 0.199 0.733 0.390 0.199 0.390 0.199 0.199 0.390 0.199 4.05  0.199 0.199 0.390
    ##  9 0.390 2.59  1.61  2.78  1.53  1.25  3.52  1.68  3.24  3.49  2.31  2.19  1.82 
    ## 10 0.390 1.61  0.733 0.199 0.390 1.44  0.199 0.199 0.199 0.199 1.75  0.733 0.199
    ## # … with 30 more rows, and 136,259 more variables: V14 <dbl>, V15 <dbl>,
    ## #   V16 <dbl>, V17 <dbl>, V18 <dbl>, V19 <dbl>, V20 <dbl>, V21 <dbl>,
    ## #   V22 <dbl>, V23 <dbl>, V24 <dbl>, V25 <dbl>, V26 <dbl>, V27 <dbl>,
    ## #   V28 <dbl>, V29 <dbl>, V30 <dbl>, V31 <dbl>, V32 <dbl>, V33 <dbl>,
    ## #   V34 <dbl>, V35 <dbl>, V36 <dbl>, V37 <dbl>, V38 <dbl>, V39 <dbl>,
    ## #   V40 <dbl>, V41 <dbl>, V42 <dbl>, V43 <dbl>, V44 <dbl>, V45 <dbl>,
    ## #   V46 <dbl>, V47 <dbl>, V48 <dbl>, V49 <dbl>, V50 <dbl>, V51 <dbl>,
    ## #   V52 <dbl>, V53 <dbl>, V54 <dbl>, V55 <dbl>, V56 <dbl>, V57 <dbl>,
    ## #   V58 <dbl>, V59 <dbl>, V60 <dbl>, V61 <dbl>, V62 <dbl>, V63 <dbl>,
    ## #   V64 <dbl>, V65 <dbl>, V66 <dbl>, V67 <dbl>, V68 <dbl>, V69 <dbl>,
    ## #   V70 <dbl>, V71 <dbl>, V72 <dbl>, V73 <dbl>, V74 <dbl>, V75 <dbl>,
    ## #   V76 <dbl>, V77 <dbl>, V78 <dbl>, V79 <dbl>, V80 <dbl>, V81 <dbl>,
    ## #   V82 <dbl>, V83 <dbl>, V84 <dbl>, V85 <dbl>, V86 <dbl>, V87 <dbl>,
    ## #   V88 <dbl>, V89 <dbl>, V90 <dbl>, V91 <dbl>, V92 <dbl>, V93 <dbl>,
    ## #   V94 <dbl>, V95 <dbl>, V96 <dbl>, V97 <dbl>, V98 <dbl>, V99 <dbl>,
    ## #   V100 <dbl>, V101 <dbl>, V102 <dbl>, V103 <dbl>, V104 <dbl>, V105 <dbl>,
    ## #   V106 <dbl>, V107 <dbl>, V108 <dbl>, V109 <dbl>, V110 <dbl>, V111 <dbl>,
    ## #   V112 <dbl>, V113 <dbl>, …

put healthy cells through pca

``` r
healthy_pcs <- 
  recipe(healthy_sampled) %>% 
  #step_center(everything()) %>% 
  #step_scale(everything()) %>% 
  step_pca(everything(), threshold = 0.95) %>% 
  prep() %>% 
  juice() %>% 
  as.matrix()

as_tibble(healthy_pcs)
```

    ## # A tibble: 40 x 10
    ##       PC01   PC02   PC03   PC04    PC05    PC06   PC07     PC08    PC09    PC10
    ##      <dbl>  <dbl>  <dbl>  <dbl>   <dbl>   <dbl>  <dbl>    <dbl>   <dbl>   <dbl>
    ##  1   -98.4   39.6   12.0   3.64   -6.63    1.24  -3.33    0.860   -1.11    7.12
    ##  2  -796.  -228.   257.  -46.2   -82.7   113.   -68.9   -21.0    -39.2    85.9 
    ##  3 -1128.  -230.   129.   36.3   184.   -221.   -74.7  -143.     -57.9   -84.6 
    ##  4  -788.  -173.    41.4 221.     94.7   245.   -96.2  -142.    -113.   -119.  
    ##  5  -737.  -400.    70.7 -29.6  -106.     30.8  -39.2   -99.1     -1.50   71.1 
    ##  6  -167.   150.   179.    2.14   61.9   -35.5  -19.7   -23.9     11.8    22.0 
    ##  7  -966.  -288.   -82.4  35.8   -41.7     8.07 -88.9     2.89    11.7   -48.7 
    ##  8  -271.   222.  -165.   81.4  -143.    143.   -13.3  -182.      -3.76    6.65
    ##  9  -685.   223.    23.4 -44.2  -151.     69.1  -58.9   121.    -126.    -86.7 
    ## 10  -289.  -103.    37.0 -44.6   -82.3     1.95 208.    -66.7     17.1   -31.4 
    ## # … with 30 more rows

create the matrix `cells` in which each row is an antigen channel and
each column is a cell from the full dataset

``` r
cells <- 
  aml_sampled %>% 
  select(any_of(healthy_colnames)) %>% 
  t()

as_tibble(cells)
```

    ## # A tibble: 40 x 22,508
    ##       V1    V2    V3    V4    V5    V6    V7    V8    V9   V10   V11   V12   V13
    ##    <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    ##  1 0.199 0.199 0.199 0.199 0.199 0.199 0.199 0.199 0.199 0.199 0.199 0.199 0.199
    ##  2 4.20  0.390 2.87  1.44  1.02  3.00  3.43  4.25  1.02  1.02  2.39  0.733 3.30 
    ##  3 4.31  4.43  3.07  4.49  0.199 4.18  3.21  4.05  2.78  4.53  0.390 1.25  3.85 
    ##  4 4.65  3.61  3.54  4.14  0.881 1.53  4.49  4.13  3.79  4.50  0.199 3.50  4.48 
    ##  5 4.05  0.199 4.44  2.52  0.390 2.42  3.85  3.72  0.199 3.44  0.199 0.199 3.90 
    ##  6 0.199 1.02  0.199 1.68  0.199 0.199 0.199 0.199 0.199 2.70  3.00  0.199 0.199
    ##  7 4.15  1.02  4.15  3.09  0.199 4.00  3.53  3.62  1.14  3.02  0.881 1.44  4.21 
    ##  8 0.733 0.390 0.733 0.199 0.199 0.199 0.390 0.390 1.02  0.199 0.199 0.390 0.390
    ##  9 2.23  0.199 2.23  0.569 0.881 2.59  1.94  1.68  2.78  0.881 3.30  2.27  1.44 
    ## 10 0.733 0.199 2.70  0.881 0.199 1.25  1.53  0.881 0.199 0.569 0.199 0.199 0.390
    ## # … with 30 more rows, and 22,495 more variables: V14 <dbl>, V15 <dbl>,
    ## #   V16 <dbl>, V17 <dbl>, V18 <dbl>, V19 <dbl>, V20 <dbl>, V21 <dbl>,
    ## #   V22 <dbl>, V23 <dbl>, V24 <dbl>, V25 <dbl>, V26 <dbl>, V27 <dbl>,
    ## #   V28 <dbl>, V29 <dbl>, V30 <dbl>, V31 <dbl>, V32 <dbl>, V33 <dbl>,
    ## #   V34 <dbl>, V35 <dbl>, V36 <dbl>, V37 <dbl>, V38 <dbl>, V39 <dbl>,
    ## #   V40 <dbl>, V41 <dbl>, V42 <dbl>, V43 <dbl>, V44 <dbl>, V45 <dbl>,
    ## #   V46 <dbl>, V47 <dbl>, V48 <dbl>, V49 <dbl>, V50 <dbl>, V51 <dbl>,
    ## #   V52 <dbl>, V53 <dbl>, V54 <dbl>, V55 <dbl>, V56 <dbl>, V57 <dbl>,
    ## #   V58 <dbl>, V59 <dbl>, V60 <dbl>, V61 <dbl>, V62 <dbl>, V63 <dbl>,
    ## #   V64 <dbl>, V65 <dbl>, V66 <dbl>, V67 <dbl>, V68 <dbl>, V69 <dbl>,
    ## #   V70 <dbl>, V71 <dbl>, V72 <dbl>, V73 <dbl>, V74 <dbl>, V75 <dbl>,
    ## #   V76 <dbl>, V77 <dbl>, V78 <dbl>, V79 <dbl>, V80 <dbl>, V81 <dbl>,
    ## #   V82 <dbl>, V83 <dbl>, V84 <dbl>, V85 <dbl>, V86 <dbl>, V87 <dbl>,
    ## #   V88 <dbl>, V89 <dbl>, V90 <dbl>, V91 <dbl>, V92 <dbl>, V93 <dbl>,
    ## #   V94 <dbl>, V95 <dbl>, V96 <dbl>, V97 <dbl>, V98 <dbl>, V99 <dbl>,
    ## #   V100 <dbl>, V101 <dbl>, V102 <dbl>, V103 <dbl>, V104 <dbl>, V105 <dbl>,
    ## #   V106 <dbl>, V107 <dbl>, V108 <dbl>, V109 <dbl>, V110 <dbl>, V111 <dbl>,
    ## #   V112 <dbl>, V113 <dbl>, …

fit linear models for sampled data

``` r
lm_result <- 
  lm(cells ~ healthy_pcs + 0)

healthy_component <- 
  lm_result$fitted.values %>% 
  t() %>% 
  as_tibble() %>% 
  rename_with(.fn = ~ str_c(., "_healthy"))

disease_component <- 
  lm_result$residuals %>% 
  t() %>% 
  as_tibble() %>% 
  rename_with(.fn = ~ str_c(., "_disease")) 
```

``` r
disease_component
```

    ## # A tibble: 22,508 x 40
    ##    CD10_disease CD11b_disease CD11c_disease CD13_disease CD14_disease
    ##           <dbl>         <dbl>         <dbl>        <dbl>        <dbl>
    ##  1     -0.0224        0.542          -0.602       0.285        0.314 
    ##  2     -0.00862      -0.134          -0.462       0.346       -0.164 
    ##  3     -0.0484       -0.354          -0.591      -0.0496       1.21  
    ##  4     -0.0574       -0.879           0.300       0.305        0.149 
    ##  5     -0.152        -0.436          -0.652       0.669        0.254 
    ##  6     -0.0621       -0.00606         0.284      -0.766       -0.113 
    ##  7     -0.0344       -0.322          -0.537       0.449        0.244 
    ##  8     -0.0470        0.0259          0.137      -0.208       -0.0677
    ##  9     -0.0259        0.104           0.381      -0.186        0.530 
    ## 10     -0.00443      -1.25           -0.385       0.268        1.08  
    ## # … with 22,498 more rows, and 35 more variables: CD16_disease <dbl>,
    ## #   CD33_disease <dbl>, CD34_disease <dbl>, CD38_disease <dbl>,
    ## #   CD41_disease <dbl>, CD45_disease <dbl>, CD45RA_disease <dbl>,
    ## #   CD47_disease <dbl>, CD49f_disease <dbl>, CD56_disease <dbl>,
    ## #   CD61_disease <dbl>, CD64_disease <dbl>, CD68_disease <dbl>,
    ## #   CD71_disease <dbl>, CD90_disease <dbl>, CD93_disease <dbl>,
    ## #   CD99_disease <dbl>, CD109_disease <dbl>, CD117_disease <dbl>,
    ## #   CD123_disease <dbl>, CD135_disease <dbl>, CCR2_disease <dbl>,
    ## #   TIM-3_disease <dbl>, HLA-DR_disease <dbl>, pAkt_disease <dbl>,
    ## #   pCreb_disease <dbl>, pErk_disease <dbl>, pS6_disease <dbl>,
    ## #   pSTAT3_disease <dbl>, pSTAT5_disease <dbl>, CEBPa_disease <dbl>,
    ## #   GATA-1_disease <dbl>, PU.1_disease <dbl>, caspase-3_disease <dbl>,
    ## #   MPO_disease <dbl>

combine new features with old features in `aml_sampled`

``` r
aml_sampled <- 
  aml_sampled %>% 
  bind_cols(healthy_component, disease_component)

aml_sampled %>% 
  colnames()
```

    ##   [1] "CD45"                        "CD61"                       
    ##   [3] "CD99"                        "127I"                       
    ##   [5] "CD45RA"                      "CD93"                       
    ##   [7] "CD3_CD19"                    "CCR2"                       
    ##   [9] "CD117"                       "CD123"                      
    ##  [11] "CD64"                        "CD90"                       
    ##  [13] "CD38"                        "CD34"                       
    ##  [15] "CEBPa"                       "pSTAT5"                     
    ##  [17] "CD11c"                       "CD13"                       
    ##  [19] "pAkt"                        "TIM-3"                      
    ##  [21] "CD56"                        "CD10"                       
    ##  [23] "PU.1"                        "CD33"                       
    ##  [25] "CD14"                        "caspase-3"                  
    ##  [27] "GATA-1"                      "pSTAT3"                     
    ##  [29] "CD41"                        "CD16"                       
    ##  [31] "CD68"                        "MPO"                        
    ##  [33] "pErk"                        "CD47"                       
    ##  [35] "CD135"                       "CD109"                      
    ##  [37] "pS6"                         "CD49f"                      
    ##  [39] "HLA-DR"                      "CD71"                       
    ##  [41] "pCreb"                       "191Ir"                      
    ##  [43] "193Ir"                       "cisplatin"                  
    ##  [45] "CD11b"                       "file_names"                 
    ##  [47] "plate"                       "patient"                    
    ##  [49] "stimulation"                 "condition"                  
    ##  [51] "cell_id"                     "cluster_flowSOM_surface"    
    ##  [53] "cluster_flowSOM_signaling"   "cluster_flowSOM_mah"        
    ##  [55] "cluster_flowSOM_cos"         "cluster_flowSOM_all_markers"
    ##  [57] "mahalanobis_HSC"             "mahalanobis_MPP"            
    ##  [59] "mahalanobis_CMP"             "mahalanobis_GMP"            
    ##  [61] "mahalanobis_MEP"             "mahalanobis_Monocyte"       
    ##  [63] "mahalanobis_DC"              "mahalanobis_Macrophage"     
    ##  [65] "mahalanobis_Thrombocyte"     "mahalanobis_cluster"        
    ##  [67] "cosine_HSC"                  "cosine_MPP"                 
    ##  [69] "cosine_CMP"                  "cosine_GMP"                 
    ##  [71] "cosine_MEP"                  "cosine_Monocyte"            
    ##  [73] "cosine_DC"                   "cosine_Macrophage"          
    ##  [75] "cosine_Thrombocyte"          "cosine_cluster"             
    ##  [77] "CD10_healthy"                "CD11b_healthy"              
    ##  [79] "CD11c_healthy"               "CD13_healthy"               
    ##  [81] "CD14_healthy"                "CD16_healthy"               
    ##  [83] "CD33_healthy"                "CD34_healthy"               
    ##  [85] "CD38_healthy"                "CD41_healthy"               
    ##  [87] "CD45_healthy"                "CD45RA_healthy"             
    ##  [89] "CD47_healthy"                "CD49f_healthy"              
    ##  [91] "CD56_healthy"                "CD61_healthy"               
    ##  [93] "CD64_healthy"                "CD68_healthy"               
    ##  [95] "CD71_healthy"                "CD90_healthy"               
    ##  [97] "CD93_healthy"                "CD99_healthy"               
    ##  [99] "CD109_healthy"               "CD117_healthy"              
    ## [101] "CD123_healthy"               "CD135_healthy"              
    ## [103] "CCR2_healthy"                "TIM-3_healthy"              
    ## [105] "HLA-DR_healthy"              "pAkt_healthy"               
    ## [107] "pCreb_healthy"               "pErk_healthy"               
    ## [109] "pS6_healthy"                 "pSTAT3_healthy"             
    ## [111] "pSTAT5_healthy"              "CEBPa_healthy"              
    ## [113] "GATA-1_healthy"              "PU.1_healthy"               
    ## [115] "caspase-3_healthy"           "MPO_healthy"                
    ## [117] "CD10_disease"                "CD11b_disease"              
    ## [119] "CD11c_disease"               "CD13_disease"               
    ## [121] "CD14_disease"                "CD16_disease"               
    ## [123] "CD33_disease"                "CD34_disease"               
    ## [125] "CD38_disease"                "CD41_disease"               
    ## [127] "CD45_disease"                "CD45RA_disease"             
    ## [129] "CD47_disease"                "CD49f_disease"              
    ## [131] "CD56_disease"                "CD61_disease"               
    ## [133] "CD64_disease"                "CD68_disease"               
    ## [135] "CD71_disease"                "CD90_disease"               
    ## [137] "CD93_disease"                "CD99_disease"               
    ## [139] "CD109_disease"               "CD117_disease"              
    ## [141] "CD123_disease"               "CD135_disease"              
    ## [143] "CCR2_disease"                "TIM-3_disease"              
    ## [145] "HLA-DR_disease"              "pAkt_disease"               
    ## [147] "pCreb_disease"               "pErk_disease"               
    ## [149] "pS6_disease"                 "pSTAT3_disease"             
    ## [151] "pSTAT5_disease"              "CEBPa_disease"              
    ## [153] "GATA-1_disease"              "PU.1_disease"               
    ## [155] "caspase-3_disease"           "MPO_disease"

compare healthy and diseased components in each channel across samples
types

``` r
univariate_plot <- function(marker_regex) { 
  hist_data <- 
    aml_sampled %>% 
    mutate(
      condition = fct_recode(condition, "cancer" = "dx", "cancer" = "rx")
    ) %>% 
    select(condition, contains(marker_regex)) %>% 
    pivot_longer(
      cols = -condition, 
      names_to = c('marker', 'component'), 
      values_to = "value", 
      names_sep = "_"
    ) 
  
  hist_medians <- 
    hist_data %>% 
    group_by(component, condition) %>% 
    summarize(value = median(value))
  
  hist_data %>% 
    ggplot(aes(x = value, fill = condition)) + 
    geom_vline(
      aes(xintercept = value, color = condition), 
      data = hist_medians
    ) + 
    geom_density(alpha = 0.4) + 
    facet_grid(rows = vars(component)) + 
    labs(subtitle = str_remove(marker_regex, "_"))
}

univariate_plot("CD45_")
```

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
univariate_plot("CD34_")
```

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
univariate_plot("CD38_")
```

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->

``` r
univariate_plot("CD34_")
```

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-4.png)<!-- -->

``` r
univariate_plot("CD56_")
```

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-5.png)<!-- -->

``` r
univariate_plot("HLA-DR_")
```

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-6.png)<!-- -->

``` r
map(str_c(healthy_colnames, "_"), univariate_plot)
```

    ## [[1]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-7.png)<!-- -->

    ## 
    ## [[2]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-8.png)<!-- -->

    ## 
    ## [[3]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-9.png)<!-- -->

    ## 
    ## [[4]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-10.png)<!-- -->

    ## 
    ## [[5]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-11.png)<!-- -->

    ## 
    ## [[6]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-12.png)<!-- -->

    ## 
    ## [[7]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-13.png)<!-- -->

    ## 
    ## [[8]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-14.png)<!-- -->

    ## 
    ## [[9]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-15.png)<!-- -->

    ## 
    ## [[10]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-16.png)<!-- -->

    ## 
    ## [[11]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-17.png)<!-- -->

    ## 
    ## [[12]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-18.png)<!-- -->

    ## 
    ## [[13]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-19.png)<!-- -->

    ## 
    ## [[14]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-20.png)<!-- -->

    ## 
    ## [[15]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-21.png)<!-- -->

    ## 
    ## [[16]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-22.png)<!-- -->

    ## 
    ## [[17]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-23.png)<!-- -->

    ## 
    ## [[18]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-24.png)<!-- -->

    ## 
    ## [[19]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-25.png)<!-- -->

    ## 
    ## [[20]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-26.png)<!-- -->

    ## 
    ## [[21]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-27.png)<!-- -->

    ## 
    ## [[22]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-28.png)<!-- -->

    ## 
    ## [[23]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-29.png)<!-- -->

    ## 
    ## [[24]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-30.png)<!-- -->

    ## 
    ## [[25]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-31.png)<!-- -->

    ## 
    ## [[26]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-32.png)<!-- -->

    ## 
    ## [[27]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-33.png)<!-- -->

    ## 
    ## [[28]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-34.png)<!-- -->

    ## 
    ## [[29]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-35.png)<!-- -->

    ## 
    ## [[30]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-36.png)<!-- -->

    ## 
    ## [[31]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-37.png)<!-- -->

    ## 
    ## [[32]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-38.png)<!-- -->

    ## 
    ## [[33]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-39.png)<!-- -->

    ## 
    ## [[34]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-40.png)<!-- -->

    ## 
    ## [[35]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-41.png)<!-- -->

    ## 
    ## [[36]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-42.png)<!-- -->

    ## 
    ## [[37]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-43.png)<!-- -->

    ## 
    ## [[38]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-44.png)<!-- -->

    ## 
    ## [[39]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-45.png)<!-- -->

    ## 
    ## [[40]]

![](disease_specific_analysis_files/figure-gfm/unnamed-chunk-11-46.png)<!-- -->

# build all-cell model with only classifier markers

# build all-cell model with only signaling markers

# build a separate linear model for each cell type
