#AML_cluster_summaries.R
#input: expression matrix of m cells by n variables, with some of the n variables denoting which clusters each cell belongs to
#output: none
#side-effects: plot cluster abundances and centroids by population, by patient, and by stimulation. Save files 
#              storing these matrices. 
#optimizations: make object-oriented, separate visualization from computing cluster matrix 
message("Summarizing clusters...")

library(plotrix)
source(file.path(CODE.DIRECTORY, "clustering", "extract_featureMatrix.R"))

#directory setup
clustering.output <- file.path(OUTPUT.DIRECTORY, "clustering_results")
load(file.path(OUTPUT.DIRECTORY, "data", "AML_matrix_clustered.RData"))

#Set up Mah.cluster as a factor variable with the appropriate labels
expression.matrix$Mah.cluster <- as.factor(expression.matrix$Mah.cluster)
levels(expression.matrix$Mah.cluster) <- CLASSIFIER.POPULATIONS

#identify cluster columns
cluster.cols <- dplyr::contains(match = "cluster", vars = colnames(expression.matrix))

#calculate total number of cells
total.cells <- nrow(expression.matrix) 

#new centroids/abundance in general
                
cluster.features <- list()
for(col in cluster.cols){
  my.name <- colnames(expression.matrix)[col]
  message("Starting ", my.name, " cluster summaries")
  feature.matrix <- extract_featureMatrix(expression.matrix = expression.matrix, 
                                          surface_markers = c(SURFACE.MARKERS, OTHER.MARKERS), 
                                          signaling_markers = SIGNALING.MARKERS, grouping_vars = NULL, 
                                          cluster_var = my.name, stim_var = NULL, central_tendency = median)
  cluster.features[[my.name]] <- feature.matrix
}

message("Summarizing by patient...")
#centroids by patient 
features.by.patient <- list()
for(col in cluster.cols){
  my.name <- colnames(expression.matrix)[col]
  message("Starting ", my.name)
  
  my.features <- extract_featureMatrix(expression.matrix = expression.matrix, surface_markers = c(SURFACE.MARKERS, OTHER.MARKERS), 
                                       signaling_markers = SIGNALING.MARKERS, grouping_vars = "patient", 
                                       cluster_var = my.name, stim_var = "stimulation", central_tendency = median)
  
  features.by.patient[[my.name]] <- my.features
}


message("Summarizing by condition...")
#features by condition
features.by.condition <- list()
for(col in cluster.cols){
  my.name <- colnames(expression.matrix)[col]
  message("Starting ", my.name)
  
  my.features <- extract_featureMatrix(expression.matrix = expression.matrix, surface_markers = c(SURFACE.MARKERS, OTHER.MARKERS), 
                                       signaling_markers = SIGNALING.MARKERS, grouping_vars = "condition", 
                                       cluster_var = my.name, stim_var = "stimulation", central_tendency = median)
  features.by.condition[[my.name]] <- my.features
}


message("Summarizing by patient and by condition...")
#features by patient/condition
features.by.pat.con <- list()

for(col in cluster.cols){
  my.name <- colnames(expression.matrix)[col]
  message("Starting ", my.name)
  
  my.features <- extract_featureMatrix(expression.matrix = expression.matrix, surface_markers = c(SURFACE.MARKERS, OTHER.MARKERS), 
                                       signaling_markers = SIGNALING.MARKERS, grouping_vars = c("patient", "condition"), 
                                       cluster_var = my.name, stim_var = "stimulation", central_tendency = median)
  
features.by.pat.con[[my.name]] <- my.features
}


#save cluster summaries
message("Saving cluster summaries!")
setwd(file.path(OUTPUT.DIRECTORY, "data"))
save(list = c("features", "features.by.condition", 
              "features.by.patient", "features.by.pat.con"), 
     file = "AML_cluster_features.RData")

#some general summary figures 



