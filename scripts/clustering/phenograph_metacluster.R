#phenograph_metacluster.R 

# Description: Performs phenograph clustering, then metaclustering on an input dataset with a specified number k of nearest 
#              neighbors in the graph construction step.
#
# Inputs:
#    expression_matrix = an [m x n] matrix of m cells and n measured signals on which the algorithm will operate. 
#    k = the number of neighbors for the phenograph algorithm to use in the initial graph construction.
#    cluster_markers = character vector of the column (marker) names to use in the clustering.
#    num_cells_mc = number of cells from each metacluster to use for monte carlo sampling and cluster metric calculations
#    make_patient_heatmap = boolean value indicating if a heatmap representing the percentage of each patient that falls in 
#    each metacluster should be plotted. Requires that expression_matrix has a column with the name "patient." Default = FALSE.
#
# Returns: A phenoMetaCluster object with the following components... 
#     k = the number of neighbors used in the clustering
#     num_clusters = an integer indicating the number of clusters that were identified by the Phenograph algorithm
#     my_clusters: an integer vector of length m indicating the cluster assignment for each cell in expression_matrix
#     my_silhouettes: a double vector of length m indicating the silhouette value for each phenograph cluster
#     patient_heatmap: a ggplot object encoding a heatmap of the distribution of each patient's cells across 
#                      each of the phenograph clusters identified in the algorithmic step. If FALSE, will not return a heatmap. 
#     heatmap_matrix: the matrix that was used to construct the patient_heatmap.
#     centroid_plot: a PCA/tSNE/UMAP plot of each cluster identified in the first stage of clustering, colored by patient. Will
#                    only plot if make_patient_heatmap is set to True. 
#     out_dir: directory in which to save the results of the analysis. If FALSE, will not save. 
#
# Optimizations: 
#
#     Currently, expression_matrix must contain a column called "patient" - change this to any categorical variable (or multiple).

source(file.path(CODE.DIRECTORY, 'helpers', 'dunn_mc.R'))
source(file.path(CODE.DIRECTORY, 'helpers', 'silhouette_mc.R'))
source(file.path(CODE.DIRECTORY, 'helpers', 'phenoPatient.R'))
library(furrr)
library(future)


phenograph_metacluster <- function(expression_matrix = NULL, k = 30, 
                               cluster_markers = colnames(expression_matrix), 
                               num_cells_mc = 1000, make_patient_heatmap = FALSE, 
                               out_dir = FALSE){
  
  patients_nested <- expression_matrix[,c(cluster_markers, "patient", "condition", "stimulation")] %>% group_by(patient, condition, stimulation) %>% nest()
  
  plan(multiprocess)
  phenograph_clusters <- future_map(.x = patients_nested$data, phenoPatient, .progress = TRUE, chosen.cols = cluster_markers, k = k) %>% 
    flatten_dbl() #NEED TO TEST TO SEE IF THIS WORKS, may not return values
  
  expression_matrix$phenograph_clusters <- phenograph_clusters
  
  #summarize cluster centroids 
  pheno_centroids <- expression_matrix %>% group_by(patient, condition, stimulation, phenograph_clusters) %>% 
    summarize_if(.predicate = is.numeric, .funs = median) %>% ungroup()
  
  pheno_numbers <- expression_matrix %>% group_by(patient, condition, stimulation, phenograph_clusters) %>% 
    dplyr::summarize(num_cells = n()) %>% ungroup()
  
  pheno_centroids <- left_join(pheno_centroids, pheno_numbers)
  
  #perform metaclustering
  metaclusters <- phenoPatient(my.df = pheno_centroids %>% ungroup(), 
                               chosen.cols = cluster_markers, k = k)
  pheno_centroids$`Phenograph Metaclusters` <- as.factor(metaclusters)
  
  centroid_means <- pheno_centroids %>% ungroup %>% dplyr::select(one_of(cluster_markers)) %>% summarize_all(median)
  
  #add phenograph metacluster column
  match_table <- pheno_centroids %>% dplyr::select(patient, condition, stimulation, phenograph_clusters, `Phenograph Metaclusters`)
  expression_matrix <- left_join(expression_matrix, match_table)

  #save num_clusters
  num_clusters <- length(unique(metaclusters))
  
  #make figures
  if(make_patient_heatmap){
    message('Plotting patient heatmap.')
    #make heatmap
    heatmap.matrix <- table(expression_matrix$`Phenograph Metaclusters`, expression_matrix$patient)
    patient.totals <- colSums(heatmap.matrix)
    heatmap.matrix2 <- sapply(X = 1:ncol(heatmap.matrix), FUN = function(x) heatmap.matrix[,x]/patient.totals[x])
    colnames(heatmap.matrix2) <- colnames(heatmap.matrix)
    patient_heatmap <- pheatmap(mat = t(heatmap.matrix2))
    
    #make pre-metacluster and metacluster pca plot (may need to play with scaling and centering)
    message('Plotting cluster PCA.')
    centroid_pca = prcomp(x = pheno_centroids %>% ungroup %>% dplyr::select(one_of(cluster_markers)), 
                          center = TRUE, scale. = FALSE)
    pca_clusters <- as_tibble(centroid_pca$x)
    pca_clusters <- cbind(pheno_centroids %>% ungroup(), pca_clusters)
    
    pca_contours <- t(t(as.matrix(expression_matrix[,cluster_markers])) - as.numeric(centroid_means))
    pca_contours <- pca_contours %*% as.matrix(centroid_pca$rotation)
    pca_contours <- as_tibble(pca_contours)
    
    centroid_plot <- ggplot(data = pca_clusters, 
                                mapping = aes(x = PC1, 
                                              y = PC2, 
                                              color = patient, 
                                              shape = `Phenograph Metaclusters`, 
                                              size = num_cells)) + 
      geom_point() + theme_bw() + geom_density2d(data = pca_contours,
                                                 mapping = aes(x = PC1, y = PC2), 
                                                 color = "black", alpha = 0.3, 
                                                 inherit.aes = FALSE)
      } else {
    patient_heatmap = NULL
    heatmap.matrix2 = NULL
    centroid_plot = NULL
  }
  
  #return silhouette values using euclidean distance
  # message('Finding average silhouette values using MC simulation')
  # my_silhouettes = silhouette_mc(expression_matrix = expression_matrix[,cluster_markers], 
  #                                clusters = expression_matrix$`Phenograph Metaclusters`, num_cells_mc = num_cells_mc, 
  #                                num_iterations = 10)
  # 
  # #return dunn values using euclidean distance
  # message('Finding average Dunn Index using MC simulation')
  # my_dunn = dunn_mc(expression_matrix = expression_matrix[,cluster_markers], 
  #                   clusters = expression_matrix$`Phenograph Metaclusters`, num_cells_mc = num_cells_mc, 
  #                   num_iterations = 10)
  
  my_silhouettes <- 0
  my_dunn <- 0
  
  
  #return final list 
  message('Writing results')
  phenoCluster_result <- list(k = k, num_clusters = num_clusters, my_clusters = expression_matrix$`Phenograph Metaclusters`, 
                              my_silhouettes = my_silhouettes, my_dunn = my_dunn, 
                              patient_heatmap = patient_heatmap, heatmap_matrix = heatmap.matrix2, 
                              centroid_plot = centroid_plot)
  class(phenoCluster_result) <- "phenoCluster"
  if(out_dir != FALSE){ 
    write_rds(x = phenoCluster_result, path = file.path(out_dir, 'metaclustering_result.rds'))
  }
  return(phenoCluster_result)
} 