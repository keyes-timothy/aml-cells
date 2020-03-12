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

library(FlowSOM)
source(file.path(CODE.DIRECTORY, 'helpers', 'dunn_mc.R'))
source(file.path(CODE.DIRECTORY, 'helpers', 'silhouette_mc.R'))

flowSOM_cluster <- function(expression_matrix = NULL, 
                               cluster_markers = colnames(expression_matrix), 
                               num_cells_mc = 1000, make_patient_heatmap = FALSE){
  
  my.FSO <- list(data = data.matrix(expression_matrix), 
                 compensate = FALSE, 
                 spillover = NULL, 
                 transform = FALSE, 
                 scale = NULL, 
                 prettyColnames = colnames(expression_matrix))
  my_SOM <- BuildSOM(fsom = my.FSO, colsToUse = which(colnames(expression_matrix) %in% CLASSIFIER.MARKERS))
  my_clusters <- my_SOM$map$mapping[,1] #flowSOM clustering must go here
  
  my_MST <- BuildMST(my_SOM, tSNE = FALSE)
  
  flowSOM_metaclusters <- MetaClustering(data = my_MST$map$codes, 
                                         method = "metaClustering_consensus", 
                                         max = 50)
  flowSOM_metaclusters_cells <- flowSOM_metaclusters[my_MST$map$mapping[,1]]
  
  
  

  #deal with cells not assigned to a cluster - take out?
  if(length(my_clusters) != nrow(expression_matrix)){ 
    cell.string <- as.character(as.numeric(1:nrow(expression_matrix)))
    bad.cells <- which(!(cell.string %in% names(my_clusters)))
    bad.cell.clusters <- rep(NA, length(bad.cells)) 
    names(bad.cell.clusters) <- as.character(bad.cells)
    my_clusters <- c(my_clusters, bad.cell.clusters)
    my_clusters <- my_clusters[order(as.numeric(names(my_clusters)))] 
  }
  my_clusters <- my_clusters %>% as.numeric()
  
  #save num_clusters
  num_clusters <- length(unique(my_clusters))
  
  #make heatmap - will only do this across patients (for now)
  if (make_patient_heatmap) { 
    heatmap.matrix <- table(my_clusters, expression_matrix$patient) 
    patient.totals <- colSums(heatmap.matrix)
    heatmap.matrix2 <- sapply(X = 1:ncol(heatmap.matrix), FUN = function(x) heatmap.matrix[,x]/patient.totals[x])
    colnames(heatmap.matrix2) <- colnames(heatmap.matrix)
    
    patient_heatmap <- pheatmap(mat = t(heatmap.matrix2), scale = "column")
  } else {
    patient_heatmap = NULL
    heatmap.matrix2 = NULL
  }
  
  my_silhouettes <- NULL
  my_dunn <- NULL
  #return silhouette values using euclidean distance
  # my_silhouettes = silhouette_mc(expression_matrix = expression_matrix[,cluster_markers], 
  #                                clusters = my_clusters, num_cells_mc = num_cells_mc, 
  #                                num_iterations = 2)
  # 
  # #return dunn values using euclidean distance
  # my_dunn = dunn_mc(expression_matrix = expression_matrix[,cluster_markers], 
  #                   clusters = my_clusters, num_cells_mc = num_cells_mc, 
  #                   num_iterations = 2)
  
  
  #return final list 
  flowCluster_result <- list(num_clusters = num_clusters, my_clusters = my_clusters, 
                              my_silhouttes = my_silhouettes, my_dunn = my_dunn, 
                              patient_heatmap = patient_heatmap, heatmap_matrix = heatmap.matrix2, 
                             my_metaclusters = flowSOM_metaclusters_cells)
  class(flowCluster_result) <- "flowCluster"
  flowCluster_result
} 

