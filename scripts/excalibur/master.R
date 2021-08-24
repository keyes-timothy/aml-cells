# running components of this project on excalibur

# run file io and basic pre-processing 
intro_input <- here::here("reports", "intro", "aml_intro_vignette.Rmd")
intro_output <- here::here("reports", "intro", "aml_intro_vignette.md")

#knitr::knit(input = intro_input, output = intro_output)

# run clustering (except phenograph)
clustering_input <- 
  here::here("reports", "clustering", "aml_clustering_vignette.Rmd")
clustering_output <- 
  here::here("reports", "clustering", "aml_clustering_vignette.md")

knitr::knit(input = clustering_input, output = clustering_output)

# run feature extraction 

# run disease-specific analysis

