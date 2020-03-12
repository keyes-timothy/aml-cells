#phenograph_cluster.R 

# Description: Performs phenograph clustering on an input dataset with a specified number k of nearest 
#              neighbors in the graph construction step. 
#
# Inputs:
#    expression_matrix = an [m x n] matrix of m cells and n measured signals on which the algorithm will operate. 
#    k = the number of neighbors for the phenograph algorithm to use in the initial graph construction.
#    cluster_markers = character vector of the column (marker) names to use in the clustering.
#    num_cells_mc = number of cells from each cluster to use for monte carlo sampling and cluster metric calculations
#    make_patient_heatmap = boolean value indicating if a heatmap representing the percentage of each patient that falls in 
#    each cluster should be plottted. Requires that expression_matrix has a column with the name "patient." Default = FALSE.
#
# Returns: A phenoCluster object with the following components... 
#     k = the number of neighbors used in the clustering
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

phenograph_cluster <- function(expression_matrix = NULL, k = 30, 
                               cluster_markers = colnames(expression_matrix), 
                               num_cells_mc = 1000, make_patient_heatmap = FALSE){ 
  result <- Rphenograph(data = (expression_matrix %>% dplyr::select(one_of(cluster_markers))), 
                        k = k)
  my_clusters <- membership(result[[2]])
  
  #deal with cells not assigned to a cluster
  if(length(my_clusters) != nrow(expression_matrix)){ 
    cell.string <- as.character(as.numeric(1:nrow(expression_matrix)))
    bad.cells <- which(!(cell.string %in% names(my_clusters)))
    bad.cell.clusters <- rep(NA, length(bad.cells)) 
    names(bad.cell.clusters) <- as.character(bad.cells)
    my_clusters <- c(my_clusters, bad.cell.clusters)
    my_clusters <- my_clusters[order(as.numeric(names(my_clusters)))] 
  }
  my_clusters <- my_clusters %>% as.numeric()
  
  #summarize clusters from first round of clustering
  
  
  #save num_clusters
  num_clusters <- length(unique(my_clusters))
  
  #make heatmap - will only do this across patients (for now)
  if(make_patient_heatmap){ 
    heatmap.matrix <- table(my_clusters, expression_matrix$patient) 
    patient.totals <- colSums(heatmap.matrix)
    heatmap.matrix2 <- sapply(X = 1:ncol(heatmap.matrix), FUN = function(x) heatmap.matrix[,x]/patient.totals[x])
    colnames(heatmap.matrix2) <- colnames(heatmap.matrix)
    
    patient_heatmap <- pheatmap(mat = t(heatmap.matrix2), scale = "column")
  } else {
    patient_heatmap = NULL
    heatmap.matrix2 = NULL
  }
  
  #return silhouette values using euclidean distance
  my_silhouettes = silhouette_mc(expression_matrix = expression_matrix[,cluster_markers], 
                                 clusters = my_clusters, num_cells_mc = num_cells_mc, 
                                 num_iterations = 2)
  
  #return dunn values using euclidean distance
  my_dunn = dunn_mc(expression_matrix = expression_matrix[,cluster_markers], 
                     clusters = my_clusters, num_cells_mc = num_cells_mc, 
                     num_iterations = 2)
  
  
  #return final list 
  phenoCluster_result <- list(k = k, num_clusters = num_clusters, my_clusters = my_clusters, 
                              my_silhouttes = my_silhouettes, my_dunn = my_dunn, 
                              patient_heatmap = patient_heatmap, heatmap_matrix = heatmap.matrix2)
  class(phenoCluster_result) <- "phenoCluster"
  phenoCluster_result
} 

