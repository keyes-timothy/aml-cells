###AML_clustering.R
###implement clustering algorithms for all AML data 


message("Setting up clustering run...")
#libraries
library(Rphenograph)
library(FlowSOM)
library(foreach)
library(doParallel)
library(purrr)
#Xshift
#hierarchical

#####   directories   #####

clustering.output <- file.path(OUTPUT.DIRECTORY, "clustering_results")

###########      phenograph      ############
num_cores <- 3 #can be changed


expression.matrix$Mah.cluster <- expression.matrix$MahID
expression.matrix$MahID <- NULL
surface.cols <- colnames(expression.matrix) %in% SURFACE.MARKERS
signaling.cols <- colnames(expression.matrix) %in% SIGNALING.MARKERS
mah.dist.cols <- dplyr::contains(match = "Mah-", vars = colnames(expression.matrix))
cluster.cols <- dplyr::contains(match = "cluster", vars = colnames(expression.matrix))
metadata.cols <- dplyr::one_of(c("stimulation", "condition", "patient"), .vars = colnames(expression.matrix))
other.cols <- !(colnames(expression.matrix) %in% c(SURFACE.MARKERS, 
                                                  colnames(expression.matrix)[mah.dist.cols],
                                                  colnames(expression.matrix)[cluster.cols], 
                                                  SIGNALING.MARKERS,
                                                  colnames(expression.matrix)[metadata.cols]))

#normalize channels - divide each cell's measurement of a given protein by the 99.5th percentile of healthy
message("Normalizing data")

max.healthies <- expression.matrix %>% ungroup() %>% dplyr::filter(condition == "Healthy") %>% 
  dplyr::select(c(which(surface.cols), which(signaling.cols), which(other.cols))) %>% 
  apply(MARGIN = 2, FUN = quantile, probs = c(0.995))

for(i in 1:length(max.healthies)){ 
  expression.matrix[,names(max.healthies)[[i]]] <- expression.matrix[,names(max.healthies)[[i]]]/max.healthies[[i]]
}


#define function to phenograph cluster a single df taken from a list of dataframes 
phenoPatient <- function(patient.name = NULL, my.df = NULL, chosen.cols = surface.cols){ 
  my.cells <- my.df %>% dplyr::filter(patient == patient.name)
  my.cells <- my.cells[,which(chosen.cols)]
  my.phenograph <- Rphenograph(data = my.cells, k = 15)
  my.clusters <- membership(my.phenograph[[2]]) %>% as.vector()
  my.cells <- cbind(my.cells, my.clusters)
  new.means <- my.cells %>% group_by(my.clusters) %>% summarize_all(mean)
  new.medians <- my.cells %>% group_by(my.clusters) %>% summarize_all(median)
  my.result = list(list(cell.assignments = my.clusters, 
                        means = as_tibble(new.means), 
                        medians = as_tibble(new.medians)))
  names(my.result) <- patient.name
  cat("\n", patient.name, " done!")
  return(my.result)
}

#run phenograph on each patient individually (in parallel)
message("Running PhenoGraph on all patient data")

patient.list <- unique(expression.matrix$patient)

my.clust <- makeCluster(num_cores, outfile = "")
registerDoParallel(my.clust)

my.phenoResult <- foreach(patient = patient.list, .combine = c, 
                          .packages = c("dplyr", "Rphenograph")) %dopar% 
  phenoPatient(patient.name = patient, my.df = expression.matrix, chosen.cols = surface.cols)

stopCluster(my.clust)

#extract relevant information from parallel computation's nested list result
my.phenoResult <- purrr::transpose(my.phenoResult)
pheno.assignments <- flatten_dbl(my.phenoResult$cell.assignments) #vector of cell assignments

#get a data matrix in which each row is a cluster identified in the single-patient phenograph runs
pheno.medians <- my.phenoResult$medians %>% enframe(name = "Patient", value = "median.tibble") %>% unnest()

#run phenograph on clusters from each patient to yield metaclusters
my.meta.phenograph <- pheno.medians %>% dplyr::select(-my.clusters, -Patient) %>% Rphenograph(k = 15)
my.metaclusters <- membership(my.meta.phenograph[[2]]) %>% as.vector()
pheno.medians <- cbind(pheno.medians, my.metaclusters)
metacluster.table <- pheno.medians %>% dplyr::select(Patient, Centroid = my.clusters, Metacluster = my.metaclusters)

#adding the phenograph metacluster assignment to expression.matrix
expression.matrix$phenograph.metacluster <- rep(0, times = nrow(expression.matrix))
for (patient in patient.list){ 
  lookup.table <- metacluster.table %>% dplyr::filter(Patient == patient) %>% dplyr::select(Centroid, Metacluster)
  row.names(lookup.table) <- lookup.table$Centroid
  dictionary <- lookup.table$Metacluster; names(dictionary) <- rownames(lookup.table)
  cell.assignments <- pheno.assignments[which(expression.matrix$patient == patient)]
  expression.matrix$phenograph.metacluster[(expression.matrix$patient == patient)] <- 
    dictionary[as.character(cell.assignments)]
} 

#phenograph cluster altogether 
#giant.phenograph <- Rphenograph(data = expression.matrix[,which(surface.cols)], k = 15)
#expression.matrix$phenograph.cluster <- membership(giant.phenograph[[2]]) %>% as.vector()

#save data structure
setwd(file.path(OUTPUT.DIRECTORY, "data"))
save(list = c("expression.matrix"), 
     file = "AML_matrix_Phenograph.RData")
setwd(CODE.DIRECTORY)

rm(lookup.table, my.meta.phenograph, metacluster.table)

#############     FlowSOM    ###############
message("FlowSOM clustering...") 
expression.matrix$patient <- as.factor(expression.matrix$patient)
expression.matrix$stimulation <- as.factor(expression.matrix$stimulation)
expression.matrix$condition <- as.factor(expression.matrix$condition)
my.FSO <- list(data = data.matrix(expression.matrix), 
               compensate = FALSE, 
               spillover = NULL, 
               transform = FALSE, 
               scale = NULL, 
               prettyColnames = colnames(expression.matrix))
my.SOM <- BuildSOM(fsom = my.FSO, colsToUse = which(surface.cols))

#are there SOM arguments that I want to pass to ... in the function above?
my.MST <- BuildMST(my.SOM, tSNE = FALSE)

# PlotStars(my.MST, backgroundValues = as.factor(my.MST$map$medianValues[,"condition"]))
# PlotStars(my.MST, view = "tSNE", backgroundValues = as.factor(my.MST$map$medianValues[,"condition"]))
# 
# PlotStars(my.MST, backgroundValues = as.factor(my.MST$map$medianValues[,"patient"]))
# PlotStars(my.MST, view = "tSNE", backgroundValues = as.factor(my.MST$map$medianValues[,"patient"]))
# 
# PlotStars(my.MST, backgroundValues = as.factor(my.MST$map$medianValues[,"phenograph.cluster"]))
# PlotStars(my.MST, view = "tSNE", backgroundValues = as.factor(my.MST$map$medianValues[,"phenograph.cluster"]))

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

expression.matrix$FlowSOM.cluster <- my.MST$map$mapping[,1]

# my.flowSOM.modes <- expression.matrix %>% dplyr::group_by(FlowSOM.cluster) %>%
#   summarize_all(Mode) 
# 
# PlotStars(my.MST, backgroundValues = as.factor(my.flowSOM.modes$condition))
# PlotStars(my.MST, view = "tSNE", backgroundValues = as.factor(my.flowSOM.modes$condition))
# 
# PlotStars(my.MST, backgroundValues = as.factor(my.flowSOM.modes$patient))
# PlotStars(my.MST, view = "tSNE", backgroundValues = as.factor(my.flowSOM.modes$patient))
# 
# PlotStars(my.MST, backgroundValues = as.factor(my.flowSOM.modes$phenograph.cluster))
# PlotStars(my.MST, view = "tSNE", backgroundValues = as.factor(my.flowSOM.modes$phenograph.cluster))
# 
# PlotStars(my.MST, backgroundValues = as.factor(my.flowSOM.modes$Mah.cluster))
# PlotStars(my.MST, view = "tSNE", backgroundValues = as.factor(my.flowSOM.modes$Mah.cluster))

flowSOM.metaclusters <- MetaClustering(data = my.MST$map$codes, 
                                       method = "metaClustering_consensus", 
                                       max = 50)
flowSOM.metaclusters.cells <- flowSOM.metaclusters[my.MST$map$mapping[,1]]

expression.matrix$FlowSOM.metacluster <- flowSOM.metaclusters.cells

setwd(file.path(OUTPUT.DIRECTORY, "data"))
save(list = c("expression.matrix"), file = "AML_matrix_FlowSOM.RData")


#################################       k-means       ###################################
message("K-means clustering...")

my.kmeans <- kmeans(x = expression.matrix[,which(surface.cols)], centers = 200, iter.max = 200,  nstart = 1)
kmeans.clusters <- my.kmeans$cluster
expression.matrix$kmeans.cluster <- kmeans.clusters

setwd(file.path(OUTPUT.DIRECTORY, "data"))
save(list = c("expression.matrix"), file = "AML_matrix_clustered.RData")

rm(my.clust, my.FSO, my.kmeans, my.MST, my.phenoResult, my.SOM, pheno.medians, 
   dictionary, cell.assignments, flowSOM.metaclusters, flowSOM.metaclusters.cells, 
   kmeans.clusters, my.metaclusters, pheno.assignments, expression.matrix)

message("Clustering analysis completed!")

